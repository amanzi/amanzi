/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "errors.hh"
#include "exceptions.hh"

#include "mfd3d.hpp"
#include "tensor.hpp"

#include "Flow_State.hpp"
#include "Darcy_PK.hpp"
#include "Matrix_MFD.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Darcy_PK::Darcy_PK(Teuchos::ParameterList& dp_list_, Teuchos::RCP<Flow_State> FS_MPC)
{
  FS = FS_MPC;
  dp_list = dp_list_;

  mesh_ = FS->get_mesh();
  dim = mesh_->space_dimension();
  for(int k=0; k<dim; k++) (*(FS->get_gravity()))[k];

  // Create the combined cell/face DoF map.
  super_map_ = create_super_map();

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

  Flow_PK::Init(FS_MPC);
}



/* ******************************************************************
* Extract information from Diffusion Problem parameter list.
****************************************************************** */
void Darcy_PK::Init(Matrix_MFD* matrix_, Matrix_MFD* preconditioner_)
{
  matrix = matrix_;
  preconditioner = preconditioner_;

  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));

  matrix = new Matrix_MFD(FS, *super_map_);

  solver = new AztecOO;
  solver->SetUserOperator(matrix);
  solver->SetPrecOperator(preconditioner);
  solver->SetAztecOption(AZ_solver, AZ_cg);

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

  // Process other fundamental structures
  K.resize(number_owned_cells);
  matrix->symbolic_assembling_global_matrices(*super_map_);

  // Preconditioner
  Teuchos::ParameterList ML_list = dp_list.sublist("ML Parameters");
  preconditioner->init_ML_preconditioner(ML_list); 
};


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Darcy_PK::process_parameter_list()
{
  Teuchos::ParameterList preconditioner_list;
  preconditioner_list = dp_list.get<Teuchos::ParameterList>("Diffusion Preconditioner");

  max_itrs = dp_list.get<int>("Max Iterations");
  err_tol = dp_list.get<double>("Error Tolerance");
}


/* ******************************************************************
*  Calculates steady-state solution assuming that abosulte permeability 
* does not depend on time.                                                    
****************************************************************** */
int Darcy_PK::advance_to_steady_state()
{
  int n = number_owned_cells + number_owned_faces;

  // define absolute permeability tensor
  std::vector<WhetStone::Tensor> K(number_owned_cells);

  double rho = *(FS->get_fluid_density());
  double mu = *(FS->get_fluid_viscosity());
  for (int c=0; c<K.size(); c++) K[c] *= rho / mu;

  // calculate and assemble elemental stifness matrices
  Epetra_Vector b((*solution).Map());

  populate_absolute_permeability_tensor(K);
  matrix->create_mfd_stiffness_matrices(K);
  matrix->assemble_global_matrices();
  matrix->apply_boundary_conditions(bc_markers, bc_values);
  matrix->computeSchurComplement();

  solver->SetRHS(&b);
  solver->SetLHS(&*solution);  // initial solution guess 

  solver->Iterate(max_itrs, err_tol);
  num_itrs = solver->NumIters();
  residual = solver->TrueResidual();

  std::cout << "Darcy solver performed " << solver->NumIters() << " iterations." << std::endl
            << "Norm of true residual = " << solver->TrueResidual() << std::endl;

  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->derive_darcy_flux(*solution, darcy_flux);  // Derive Darcy fluxes on faces.
  return 0;
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Darcy_PK::populate_absolute_permeability_tensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& permeability = FS->ref_absolute_permeability();

  for (int c=cmin; c<=cmax; c++) {
    K[c].init(dim, 1);
    K[c](0, 0) = permeability[c];
  }
}


/* ******************************************************************
*  Printing information about Flow status.                                                     
****************************************************************** */
void Darcy_PK::print_statistics() const
{
  if (!MyPID && verbosity_level > 0) {
    cout << "Flow PK:" << endl;
    cout << "    Execution mode = " << (standalone_mode ? "standalone" : "MPC") << endl;
    cout << "    Verbosity level = " << verbosity_level << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


/*
int Darcy_PK::Init()
{
  matvec_ = new DarcyMatvec(this);

  rho_ = list.get<double>("fluid density");
  mu_ = list.get<double>("fluid viscosity");
  gravity_ = list.get<double>("gravity");
  g_[0] = 0.0; 
  g_[1] = 0.0; 
  g_[2] = -gravity_;
  
  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(list.sublist("boundary conditions",true));
  FlowBCFactory bc_factory(mesh, bc_list);
  bc_press_ = bc_factory.CreatePressure();
  bc_head_  = bc_factory.CreateStaticHead(0.0, rho_, gravity_);
  bc_flux_  = bc_factory.CreateMassFlux();
  validate_boundary_conditions();
  
  // Evaluate the BC at dummy time 0.
  //TODO: introduce some method to have a time-parametrized, steady Darcy problem.
  bc_press_->Compute(0.0);
  bc_head_->Compute(0.0);
  bc_flux_->Compute(0.0);

  // Create the diffusion matrix (structure only, no values)
  D_ = Teuchos::rcp<DiffusionMatrix>(create_diff_matrix_(mesh));

  // Create the preconditioner (structure only, no values)
  Teuchos::ParameterList diffprecon_plist = list.sublist("Diffusion Preconditioner");
  precon_ = new DiffusionPrecon(D_, diffprecon_plist, Map());

  // Create the RHS vector.
  rhs_ = new Epetra_Vector(Map());

  // DEFINE DEFAULTS FOR PROBLEM PARAMETERS
  k_.resize(CellMap(true).NumMyElements(), 1.0);
}


DiffusionMatrix* DarcyProblem::create_diff_matrix_(const Teuchos::RCP<AmanziMesh::Mesh> &mesh) const
{
  // Generate the list of all Dirichlet-type faces.
  // The provided list should include USED BC faces.
  std::vector<int> dir_faces;
  for (BoundaryFunction::Iterator i = bc_press_->begin(); i != bc_press_->end(); ++i)
    dir_faces.push_back(i->first);
  for (BoundaryFunction::Iterator i = bc_head_->begin();  i != bc_head_->end();  ++i)
    dir_faces.push_back(i->first);
  return new DiffusionMatrix(mesh, dir_faces);
}


void DarcyProblem::set_absolute_permeability(const Epetra_Vector &k)
{
  // should verify k.Map() is CellMap()
  // should verify k values are all > 0
  Epetra_Vector k_ovl(mesh_->cell_map(true));
  Epetra_Import importer(mesh_->cell_map(true), mesh_->cell_map(false));

  k_ovl.Import(k, importer, Insert);

  for (int i=0; i<k_.size(); ++i) k_[i] = k_ovl[i];
}


void DarcyProblem::ComputeF(const Epetra_Vector &X, Epetra_Vector &F)
{
  // The cell and face-based DoF are packed together into the X and F Epetra
  // vectors: cell-based DoF in the first part, followed by the face-based DoF.
  // In addition, only the owned DoF belong to the vectors.

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
  //bc_press_->Compute(time);
  std::map<int,double> Fpress;
  for (BoundaryFunction::Iterator bc = bc_press_->begin(); bc != bc_press_->end(); ++bc) {
    Fpress[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }
  //bc_head_->Compute(time);
  std::map<int,double> Fhead;
  for (BoundaryFunction::Iterator bc = bc_head_->begin(); bc != bc_head_->end(); ++bc) {
    Fhead[bc->first] = Pface[bc->first] - bc->second;
    Pface[bc->first] = bc->second;
  }

  // Create cell and face result vectors that include ghosts.
  Epetra_Vector Fcell(CellMap(true));
  Epetra_Vector Fface(FaceMap(true));

  int cface[6];
  double aux1[6], aux2[6], gflux[6];

  Fface.PutScalar(0.0);
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    double K = (rho_ * k_[j] / mu_);
    MD[j].diff_op(K, Pcell[j], aux1, Fcell[j], aux2);
    // Gravity contribution
    MD[j].GravityFlux(g_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] -= rho_ * K * gflux[k];
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) Fface[cface[k]] += aux2[k];
  }

  // Dirichlet-type condition residuals; overwrite with the pre-computed values.
  std::map<int,double>::const_iterator i;
  for (i = Fpress.begin(); i != Fpress.end(); ++i) Fface[i->first] = i->second;
  for (i = Fhead.begin();  i != Fhead.end();  ++i) Fface[i->first] = i->second;
  
  // Mass flux BC contribution.
  //bc_flux_->Compute(time);
  for (BoundaryFunction::Iterator bc = bc_flux_->begin(); bc != bc_flux_->end(); ++bc)
    Fface[bc->first] += bc->second * md_->face_area_[bc->first];
  
  // Copy owned part of result into the output vectors.
  for (int j = 0; j < Fcell_own.MyLength(); ++j) Fcell_own[j] = Fcell[j];
  for (int j = 0; j < Fface_own.MyLength(); ++j) Fface_own[j] = Fface[j];

  delete &Pcell_own, &Pface_own, &Fcell_own, &Fface_own;
}


Epetra_Vector* DarcyProblem::CreateCellView(const Epetra_Vector &X) const
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  return new Epetra_Vector(View, CellMap(), data);
}


Epetra_Vector* DarcyProblem::CreateFaceView(const Epetra_Vector &X) const
{
  // should verify that X.Map() is the same as Map()
  double *data;
  X.ExtractView(&data);
  int ncell = CellMap().NumMyElements();
  return new Epetra_Vector(View, FaceMap(), data+ncell);
}


void DarcyProblem::DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const
{
  // should verify X.Map() is Map()
  // should verify Q.Map() is CellMap(false)
  // should verify Q.NumVectors() is 3

  // Create views into the cell and face segments of X and F
  Epetra_Vector &Pcell = *CreateCellView(X);
  Epetra_Vector &Pface_own = *CreateFaceView(X);

  // Create face vectors that include ghosts.
  Epetra_Vector Pface(FaceMap(true));

  // Populate the face pressure vector from the input.
  Pface.Import(Pface_own, *face_importer_, Insert);

  int cface[6];
  double aux1[6], aux2[6], aux3[3], gflux[6], dummy;

  for (int j = 0; j < CellMap(false).NumMyElements(); ++j) {

  }
  for (int j = 0; j < Pcell.MyLength(); ++j) {
    mesh_->cell_to_faces((unsigned int) j, (unsigned int*) cface, (unsigned int*) cface+6);
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    double K = (k_[j] / mu_);
    MD[j].diff_op(K, Pcell[j], aux1, dummy, aux2);
    MD[j].GravityFlux(g_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] = rho_ * K * gflux[k] - aux2[k];
    MD[j].CellFluxVector(aux2, aux3);
    Q[0][j] = aux3[0];
    Q[1][j] = aux3[1];
    Q[2][j] = aux3[2];
  }

  delete &Pcell, &Pface_own;
}


void DarcyProblem::DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l2_error) const
{
  /// should verify P.Map() is Map()
  /// should verify F.Map() is FaceMap()

  int fdirs[6];
  unsigned int cface[6];
  double aux1[6], aux2[6], gflux[6], dummy;

  // Create a view into the cell pressure segment of P.
  Epetra_Vector &Pcell_own = *CreateCellView(P);

  // Create a copy of the face pressure segment of P that includes ghosts.
  Epetra_Vector &Pface_own = *CreateFaceView(P);
  Epetra_Vector Pface(FaceMap(true));
  Pface.Import(Pface_own, *face_importer_, Insert);

  // Create face flux and face count vectors that include ghosts.
  Epetra_Vector Fface(FaceMap(true)); // fills with 0
  Epetra_Vector error(FaceMap(true)); // fills with 0
  Epetra_IntVector count(FaceMap(true)); // fills with 0

  // Process-local assembly of the face mass fluxes.
  for (unsigned int j = 0; j < Pcell_own.MyLength(); ++j) {
    // Get the list of process-local face indices for this cell.
    mesh_->cell_to_faces(j, cface, cface+6);
    // Gather the local face pressures int AUX1.
    for (int k = 0; k < 6; ++k) aux1[k] = Pface[cface[k]];
    // Compute the local value of the diffusion operator.
    double K = (rho_ * k_[j] / mu_);
    MD[j].diff_op(K, Pcell_own[j], aux1, dummy, aux2);
    // Gravity contribution
    MD[j].GravityFlux(g_, gflux);
    for (int k = 0; k < 6; ++k) aux2[k] = rho_ * K * gflux[k] - aux2[k];
    mesh_->cell_to_face_dirs(j, fdirs, fdirs+6);
    // Scatter the local face result into FFACE.
    for (int k = 0; k < 6; ++k) {
      Fface[cface[k]] += fdirs[k] * aux2[k];
      error[cface[k]] += aux2[k]; // sums to the flux discrepancy
      count[cface[k]]++;
    }
  }

  // Global assembly of face mass fluxes into the return vector.
  F.Export(Fface, *face_importer_, Add);

  // Create an owned face count vector that overlays the count vector with ghosts.
  int *count_data;
  count.ExtractView(&count_data);
  Epetra_IntVector count_own(View, FaceMap(false), count_data);

  // Final global assembly of the face count vector.
  count_own.Export(count, *face_importer_, Add);

  // Correct the double counting of fluxes on interior faces and convert
  // mass flux to Darcy flux by dividing by the constant fluid density.
  for (int j = 0; j < F.MyLength(); ++j)
    F[j] = F[j] / (rho_ * count[j]);

  // Create an owned face error vector that overlays the error vector with ghosts.
  double *error_data;
  error.ExtractView(&error_data);
  Epetra_Vector error_own(View, FaceMap(false), error_data);

  // Final global assembly of the flux discrepancy vector.
  error_own.Export(error, *face_importer_, Add);

  // Set the flux discrepancy error to 0 on boundary faces where there was only one flux computed.
  for (int j = 0; j < F.MyLength(); ++j)
    if (count[j] == 1) error_own[j] = 0.0;

  // Compute the norm of the flux discrepancy.
  error_own.Norm2(&l2_error);

  delete &Pcell_own, &Pface_own;
}
*/


