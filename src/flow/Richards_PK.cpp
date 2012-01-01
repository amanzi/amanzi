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

#include "Flow_BC_Factory.hpp"
#include "boundary-function.hh"
#include "Richards_PK.hpp"
#include "Interface_BDF2.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::RCP<Teuchos::ParameterList> rp_list_, Teuchos::RCP<Flow_State> FS_MPC)
{
  Flow_PK::Init(FS_MPC);

  FS = FS_MPC;
  rp_list = rp_list_;

  mesh_ = FS->get_mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = create_super_map();

  // Other fundamental physical quantaties
  rho = *(FS->get_fluid_density());
  mu = *(FS->get_fluid_viscosity()); 
  gravity.init(dim);
  for (int k=0; k<dim; k++) gravity[k] = (*(FS->get_gravity()))[k];

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

  solver = new AztecOO;
  solver->SetUserOperator(matrix);
  solver->SetPrecOperator(preconditioner);
  solver->SetAztecOption(AZ_solver, AZ_cg);

  // Get some solver parameters from the flow parameter list.
  processParameterList();

  // Process boundary data
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers.resize(nfaces, FLOW_BC_FACE_NULL);
  bc_values.resize(nfaces, 0.0);

  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;

  bc_pressure->Compute(time);
  updateBoundaryConditions(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

  // Process other fundamental structures
  K.resize(number_owned_cells);
  matrix->symbolicAssembleGlobalMatrices(*super_map_);

  // Preconditioner
  Teuchos::ParameterList ML_list = rp_list->sublist("ML Parameters");
  preconditioner->init_ML_preconditioner(ML_list); 

  // Create the Richards model evaluator and time integrator
  Teuchos::ParameterList& rme_list = rp_list->sublist("Richards model evaluator");
  solver_BDF = new Interface_BDF2(this, rme_list);  

  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(rp_list->sublist("Time integrator")));
  time_stepper = new BDF2::Dae(*solver_BDF, *super_map_);
  time_stepper->setParameterList(bdf2_list);

  // Allocate data for relative permeability
  Krel_cells = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(true)));
  if (upwind_Krel) Krel_faces = Teuchos::rcp(new Epetra_Vector(mesh_->face_map(true)));
};


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advanceSteadyState_Picard()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;

  //solver->SetAztecOption(AZ_conv, AZ_rhs);
  int itrs = 0;
  double L2norm, L2error = 1.0, L2error_prev = 1.0;
  double relaxation = 0.8;

  while (L2error > err_tol_sss && itrs < max_itrs_sss) {
    // work-around limited support for tensors
    setAbsolutePermeabilityTensor(K);
    if (upwind_Krel) calculateRelativePermeabilityUpwindGravity(*solution_cells);
    else calculateRelativePermeability(*solution_cells);

    for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;
    //for (int c=0; c<K.size(); c++) K[c] *= rho / mu;

    matrix->createMFDstiffnessMatrices(K);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    rhs = matrix->get_rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, err_tol);
    num_itrs = solver->NumIters();
    residual = solver->TrueResidual();

    L2norm = L2error = 0.0;
    for (int c=0; c<number_owned_cells; c++) {
      solution_new[c] = relaxation * solution_old[c] + (1.0 - relaxation) * solution_new[c];
      L2norm += pow(solution_new[c], 2.0);
      L2error += pow(solution_old[c] - solution_new[c], 2.0);
      solution_old[c] = solution_new[c];
    }
    L2error = sqrt(L2error / L2norm);

    if (L2error > L2error_prev) relaxation = std::min(0.95, relaxation * 1.5); 
    else relaxation = std::max(0.05, relaxation * 0.9);  
    cout << itrs << " ERR=" << L2error << " relax=" << relaxation << endl;

    L2error_prev = L2error;
    itrs++;
  }    
//for (int c=0; c<K.size(); c+=4) cout << solution_new[c] << endl; 
  
  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->createMFDstiffnessMatrices(K);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyFlux(*solution, *face_importer_, darcy_flux);
  //addGravityFluxes(darcy_flux);

  return 0;
}


/* ******************************************************************* 
 * . 
 ****************************************************************** */
int Richards_PK::advance(double dT) 
{
  // work-around limited support for tensors
  std::vector<WhetStone::Tensor> K(number_owned_cells);
  setAbsolutePermeabilityTensor(K);
  calculateRelativePermeability(*solution_cells);

  for (int c=0; c<K.size(); c++) K[c] *= rho / mu;

  matrix->createMFDstiffnessMatrices(K);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();
  matrix->computeSchurComplement(bc_markers, bc_values);
  matrix->update_ML_preconditioner();

  int itrs;
  double T0, dTnext;
  if (itrs == 0) {  // initialization
    Epetra_Vector pdot(*super_map_);
    computePDot(T0, *solution, pdot);
    time_stepper->set_initial_state(T0, *solution, pdot);

    int errc;
    double dT0 = 1e-3;
    solver_BDF->update_precon(T0, *solution, dT0, errc);
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
* .                                               
****************************************************************** */
void Richards_PK::computeFunctionalPTerm(const Epetra_Vector &X, Epetra_Vector &F, double time)
{
  // Create views into the cell and face segments of X and F
  /*
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

  // Compute the functional  
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
  */
}


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
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
void Richards_PK::calculateRelativePermeability(const Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();

    AmanziMesh::Set_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) (*Krel_cells)[*i] = WRM[mb]->k_relative(p[*i]);
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p)
{
  calculateRelativePermeability(p);  // populates cell-based permeabilities
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
void Richards_PK::derivedSdP(const Epetra_Vector& p, Epetra_Vector& ds)
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
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) s[*i] = WRM[mb]->saturation(p[*i]);
  }
}


/* ******************************************************************
* .                                               
****************************************************************** */
void Richards_PK::derivePressureFromSaturation(double s, Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) p[*i] = WRM[mb]->pressure(s);
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
* Routine updates elemental discretization matrices and must be 
* called before applying boundary conditions and global assembling.                                             
****************************************************************** */
void Richards_PK::addGravityFluxes_MFD(Matrix_MFD* matrix)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  std::vector<double> gravity_flux;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    mesh_->cell_get_face_dirs(c, &dirs);
    int nfaces = faces.size();

    calculateGravityFluxes(c, K[c], gravity_flux);
    
    Epetra_SerialDenseVector& Ff = matrix->get_Ff_cells()[c];
    for (int n=0; n<nfaces; n++) Ff[n] += gravity_flux[n] * dirs[n];
  }
}


/* ******************************************************************
* Adds .                                               
****************************************************************** */
void Richards_PK::addTimeDerivative(Epetra_Vector& pressure_cells, Matrix_MFD* matrix)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  derivedSdP(pressure_cells, dSdP);
 
  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    Acc_cells[c] += dSdP[c] * phi[c] * volume / dT;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

