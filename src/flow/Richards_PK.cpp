/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>

#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"

#include "Richards_PK.hpp"
#include "Flow_BC_Factory.hpp"
#include "boundary-function.hh"
#include "vanGenuchtenModel.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& rp_list_, Teuchos::RCP<Flow_State> FS_MPC)
{
  FS = FS_MPC;
  rp_list = rp_list_;

  mesh_ = FS->get_mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = create_super_map();

  // Other fundamental physical quantaties
  double rho = *(FS->get_fluid_density());
  double mu = *(FS->get_fluid_viscosity()); 
  for(int k=0; k<dim; k++) (*(FS->get_gravity()))[k];

#ifdef HAVE_MPI
  const  Epetra_Comm & comm = mesh_->cell_map(false).Comm(); 
  MyPID = comm.MyPID();

  const Epetra_Map& source_cmap = mesh_->cell_map(false);
  const Epetra_Map& target_cmap = mesh_->cell_map(true);

  cell_importer_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));

  const Epetra_Map& source_fmap = mesh_->face_map(false);
  const Epetra_Map& target_fmap = mesh_->face_map(true);

  face_importer_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
#endif

  // miscalleneous
  upwind_Krel = true;

  Flow_PK::Init(FS_MPC);
}


/* ******************************************************************
* Extract information from Richards Problem parameter list.
****************************************************************** */
void Richards_PK::Init(Matrix_MFD* matrix_, Matrix_MFD* preconditioner_)
{
  if (matrix_ == NULL) matrix = new Matrix_MFD(FS, *super_map_);
  else matrix = matrix_;

  if (preconditioner_ == NULL) preconditioner = matrix;
  else preconditioner = preconditioner_;

  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));

  // Create the solution vectors.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->create_cell_view(*solution));
  solution_faces = Teuchos::rcp(FS->create_face_view(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));

  // Get some solver parameters from the flow parameter list.
  process_parameter_list();

  // Process boundary data
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  std::vector<int> bc_markers(FLOW_BC_FACE_NULL, ncells);
  std::vector<double> bc_values(0.0, ncells);

  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;

  bc_pressure->Compute(time);
  update_data_boundary_faces(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

  // Create the Richards model evaluator and time integrator
  Teuchos::ParameterList& rme_list = rp_list.sublist("Richards model evaluator");
  solver = new Interface_BDF2();  

  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(rp_list.sublist("Time integrator")));
  time_stepper = new BDF2::Dae(*solver, *super_map_);
  time_stepper->setParameterList(bdf2_list);

  // Allocate data for relative permeability
  Krel_cells = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(true)));
  if (upwind_Krel) Krel_faces = Teuchos::rcp(new Epetra_Vector(mesh_->face_map(true)));
};


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Richards_PK::process_parameter_list()
{
  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = rp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  max_itrs = rp_list.get<int>("Max Iterations");
  err_tol = rp_list.get<double>("Error Tolerance");

  upwind_Krel = rp_list.get<bool>("Upwind relative permeability", true);
  atm_pressure = rp_list.get<double>("Atmospheric pressure");

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(rp_list.sublist("boundary conditions", true));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.CreatePressure();
  bc_head = bc_factory.CreateStaticHead(0.0, rho, gravity[dim - 1]);
  bc_flux = bc_factory.CreateMassFlux();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);  

  double time = (standalone_mode) ? T_internal : T_physical;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);

  // Create water retention models
  if (!rp_list.isSublist("Water retention models")) {
    Errors::Message m("There is no Water retention models list");
    Exceptions::amanzi_throw(m);
  }
  Teuchos::ParameterList &vG_list = rp_list.sublist("Water retention models");

  int nblocks = 0;  // Find out how many WRM entries there are.
  for (Teuchos::ParameterList::ConstIterator i = vG_list.begin(); i != vG_list.end(); i++) {
    if (vG_list.isSublist(vG_list.name(i))) {
      nblocks++;
    } else {
      Errors::Message msg("Water retention models sublist contains an entry that is not a sublist.");
      Exceptions::amanzi_throw(msg);
    }
  }

  WRM.resize(nblocks);
  
  int iblock = 0;
  for (Teuchos::ParameterList::ConstIterator i = vG_list.begin(); i != vG_list.end(); i++) {
    if (vG_list.isSublist(vG_list.name(i))) {
      Teuchos::ParameterList &wrm_list = vG_list.sublist(vG_list.name(i));

      if ( wrm_list.get<string>("Water retention model") == "van Genuchten") {
        std::string region = wrm_list.get<std::string>("Region");  // associated mesh block

        double vG_m = wrm_list.get<double>("van Genuchten m");
        double vG_alpha = wrm_list.get<double>("van Genuchten alpha");
        double vG_sr = wrm_list.get<double>("van Genuchten residual saturation");
	      
        WRM[iblock] = Teuchos::rcp(new vanGenuchtenModel(region, vG_m, vG_alpha, vG_sr, atm_pressure));
      }
      iblock++;
    }
  }
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advance_to_steady_state()
{
  double T0 = (standalone_mode) ? T_internal : T_physical;
  double dTnext, T1, Tlast = T0;

  int itrs = 0;
  while (T1 < Tlast) {
    // work-around limited support for tensors
    std::vector<WhetStone::Tensor> K(number_owned_cells);
    populateAbsolutePermeabilityTensor(K);
    if (upwind_Krel) relativePermeabilityUpwindGravity(*solution_cells);
    else relativePermeability(*solution_cells);

    for (int c=0; c<K.size(); c++) K[c] *= rho / mu;

    matrix->createMFDstiffnessMatrices(K);
    matrix->assembleGlobalMatrices(*rhs);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->computeSchurComplement();

    if (itrs == 0) {  // initialization
      Epetra_Vector pdot(*super_map_);
      computePDot(T0, *solution, pdot);
      time_stepper->set_initial_state(T0, *solution, pdot);

      int errc;
      double dT0 = 1e-3;
      solver->update_precon(T0, *solution, dT0, errc);
    }

    time_stepper->bdf2_step(dT, 0.0, *solution, dTnext);
    time_stepper->commit_solution(dT, *solution);
    time_stepper->write_bdf2_stepping_statistics();

    dT = dTnext;
    Tlast = time_stepper->most_recent_time();
    itrs++;
  }    
  
  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->deriveDarcyFlux(*solution, *rhs, *face_importer_, darcy_flux);
}


/* ******************************************************************* 
 * . 
 ****************************************************************** */
int Richards_PK::advance(double dT) 
{
  // work-around limited support for tensors
  std::vector<WhetStone::Tensor> K(number_owned_cells);
  populateAbsolutePermeabilityTensor(K);
  relativePermeability(*solution_cells);

  for (int c=0; c<K.size(); c++) K[c] *= rho / mu;

  matrix->createMFDstiffnessMatrices(K);
  matrix->assembleGlobalMatrices(*rhs);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->computeSchurComplement();

  int itrs;
  double T0, dTnext;
  if (itrs == 0) {  // initialization
    Epetra_Vector pdot(*super_map_);
    computePDot(T0, *solution, pdot);
    time_stepper->set_initial_state(T0, *solution, pdot);

    int errc;
    double dT0 = 1e-3;
    solver->update_precon(T0, *solution, dT0, errc);
  }

  time_stepper->bdf2_step(dT, 0.0, *solution, dTnext);
  time_stepper->commit_solution(dT, *solution);
  time_stepper->write_bdf2_stepping_statistics();
}


/* ******************************************************************
* Estimate dp/dt from the pressure equations.                                               
****************************************************************** */
void Richards_PK::computePDot(const double T, const Epetra_Vector& p, Epetra_Vector &pdot)
{
  computeFunctionalPTerm(p, pdot);

  // zero out the face part, why? (lipnikov@lanl.gov)
  Epetra_Vector *pdot_face = FS->create_face_view(pdot);
  pdot_face->PutScalar(0.0);
}


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::populateAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& permeability = FS->ref_absolute_permeability();

  for (int c=cmin; c<=cmax; c++) {
    K[c].init(dim, 1);
    K[c](0, 0) = permeability[c];
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::relativePermeability(const Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::USED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::USED, &block);

    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) (*Krel_cells)[*i] = WRM[mb]->k_relative(p[*i]);
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::relativePermeabilityUpwindGravity(const Epetra_Vector& p)
{
  relativePermeability(p);  // populates cell-based permeabilities
  FS->distribute_cell_vector(*Krel_cells);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c=0; c<number_owned_cells; c++) {
    mesh_->cell_get_face_dirs(c, &dirs);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n<nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f); 
      if ((normal * gravity) * dirs[n] > 0.0) (*Krel_faces)[f] = (*Krel_cells)[c]; 
    }
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::dSofP(const Epetra_Vector& p, Epetra_Vector& ds)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) ds[*i] = WRM[mb]->d_saturation(p[*i]);
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::deriveVanGenuchtenSaturation(const Epetra_Vector& p, Epetra_Vector& s)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, manziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) S[*i] = WRM[mb]->saturation(p[*i]);
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::setInitialPressureFromSaturation(double s, Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); ++) p[*i] = WRM[mb]->pressure(s);
  } 
 
  /*
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) p[*i] = WRM[mb]->pressure(s);
  }
  */
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::ComputePrecon(const Epetra_Vector &P, const double h)
{
  std::vector<double> K(k_);

  Epetra_Vector &Pcell_own = *CreateCellView(P);
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);

  Epetra_Vector &Pface_own = *CreateFaceView(P);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);
  
  if (upwind_k_rel_) {
    Epetra_Vector K_upwind(Pface.Map());
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int j = 0; j < K.size(); ++j) K[j] = rho_ * k_[j] / mu_;
    D_->Compute(K, K_upwind);
  } else {
    Epetra_Vector k_rel(Pcell.Map());
    ComputeRelPerm(Pcell, k_rel);
    for (int j = 0; j < K.size(); ++j) K[j] = rho_ * k_[j] * k_rel[j] / mu_;
    D_->Compute(K);
  }

  // add the time derivative to the diagonal
  Epetra_Vector celldiag(CellMap(false));
  dSofP(Pcell_own, celldiag);
 
  // get the porosity 
  const Epetra_Vector& phi = FS->get_porosity();

  celldiag.Multiply(rho_, celldiag, phi, 0.0);
   celldiag.Multiply(1.0/h, celldiag, *cell_volumes, 0.0);
  
  D_->add_to_celldiag(celldiag);

  // Compute the face Schur complement of the diffusion matrix.
  D_->ComputeFaceSchur();

  // Compute the preconditioner from the newly computed diffusion matrix and Schur complement.
  precon_->Compute();

  delete &Pcell_own, &Pface_own;
}

}  // namespace AmanziFlow
}  // namespace Amanzi


/*
void RichardsProblem::ComputeF(const Epetra_Vector &X, Epetra_Vector &F, double time)
{
  // Create views into the cell and face segments of X and F
  Epetra_Vector &Pcell_own = *CreateCellView(X);
  Epetra_Vector &Pface_own = *CreateFaceView(X);

  Epetra_Vector &Fcell_own = *CreateCellView(F);
  Epetra_Vector &Fface_own = *CreateFaceView(F);

  // Create input cell and face pressure vectors that include ghosts.
  Epetra_Vector Pcell(CellMap(true));
  Pcell.Import(Pcell_own, *cell_importer_, Insert);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Precompute the Dirichlet-type BC residuals for later use.
  // Impose the Dirichlet boundary data on the face pressure vector.
  bc_press_->Compute(time);
  std::map<int,double> Fpress;
  for (BoundaryFunction::Iterator bc = bc_press_->begin(); bc != bc_press_->end(); ++bc) {
    Fpress[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }
  bc_head_->Compute(time);
  std::map<int,double> Fhead;
  for (BoundaryFunction::Iterator bc = bc_head_->begin(); bc != bc_head_->end(); ++bc) {
    Fhead[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }

  // Computate the functional  
  Epetra_Vector K(CellMap(true)), K_upwind(FaceMap(true));
  if (upwind_k_rel_) {
    ComputeUpwindRelPerm(Pcell, Pface, K_upwind);
    for (int c=0; c<K.MyLength(); ++c) K[c] = rho_ * k_[c] / mu_;
  } else {
    ComputeRelPerm(Pcell, K);
    for (int c=0; c<K.MyLength(); ++c) K[c] = rho_ * k_[c] * K[c] / mu_;
  }

  // Create cell and face result vectors that include ghosts.
  Epetra_Vector Fcell(CellMap(true));
  Epetra_Vector Fface(FaceMap(true));

  int cface[6];
  double aux1[6], aux2[6], aux3[6], gflux[6];

  Fface.PutScalar(0.0);
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    if (upwind_k_rel_) {
      for (int k=0; k<6; ++k) aux3[k] = K_upwind[cface[k]];
      MD[j].diff_op(K[j], aux3, Pcell[j], aux1, Fcell[j], aux2);
      // Gravity contribution
      MD[j].GravityFlux(gvec_, gflux);
      for (int k = 0; k < 6; ++k) gflux[k] *= rho_ * K[j] * aux3[k];
      for (int k = 0; k < 6; ++k) {
        aux2[k] -= gflux[k];
        Fcell[j] += gflux[k];
      }
    } else {
      MD[j].diff_op(K[j], Pcell[j], aux1, Fcell[j], aux2);
      // Gravity contribution
      MD[j].GravityFlux(gvec_, gflux);
      for (int k = 0; k < 6; ++k) aux2[k] -= rho_ * K[j] * gflux[k];
    }
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) Fface[cface[k]] += aux2[k];
  }

  // Dirichlet-type condition residuals; overwrite with the pre-computed values.
  std::map<int,double>::const_iterator i;
  for (i = Fpress.begin(); i != Fpress.end(); ++i) Fface[i->first] = i->second;
  for (i = Fhead.begin();  i != Fhead.end();  ++i) Fface[i->first] = i->second;
  
  // Mass flux BC contribution.
  bc_flux_->Compute(time);
  for (BoundaryFunction::Iterator bc = bc_flux_->begin(); bc != bc_flux_->end(); ++bc)
    Fface[bc->first] += bc->second * md_->face_area_[bc->first];
  
  // Copy owned part of result into the output vectors.
  for (int j = 0; j < Fcell_own.MyLength(); ++j) Fcell_own[j] = Fcell[j];
  for (int j = 0; j < Fface_own.MyLength(); ++j) Fface_own[j] = Fface[j];

  delete &Pcell_own, &Pface_own, &Fcell_own, &Fface_own;
}
*/
