/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

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
* each variable initialization
****************************************************************** */
Darcy_PK::Darcy_PK(Teuchos::ParameterList& global_list, Teuchos::RCP<Flow_State> FS_MPC)
{
  Flow_PK::Init(FS_MPC);  // sets up default parameters
  FS = FS_MPC;

  // extract important sublists
  Teuchos::ParameterList flow_list;
  if (global_list.isSublist("Flow")) {
    flow_list = global_list.sublist("Flow");
  } else {
    Errors::Message msg("Darcy PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Darcy Problem")) {
    dp_list_ = flow_list.sublist("Darcy Problem");
  } else {
    Errors::Message msg("Darcy PK: input parameter list does not have <Darcy Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list.isSublist("Preconditioners")) {
    preconditioner_list_ = global_list.sublist("Preconditioners");
  } else {
    Errors::Message msg("Darcy PK: input parameter list does not have <Preconditioners> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list.isSublist("Solvers")) {
    solver_list_ = global_list.sublist("Solvers");
  } else {
    Errors::Message msg("Darcy PK: input parameter list does not have <Solvers> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  mesh_ = FS->mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = CreateSuperMap();

  // Other fundamental physical quantities
  rho_ = *(FS->fluid_density());
  mu_ = *(FS->fluid_viscosity());
  gravity_.init(dim);
  for (int k = 0; k < dim; k++) gravity_[k] = (*(FS->gravity()))[k];

#ifdef HAVE_MPI
  const  Epetra_Comm& comm = mesh_->cell_map(false).Comm();
  MyPID = comm.MyPID();

  const Epetra_Map& source_cmap = mesh_->cell_map(false);
  const Epetra_Map& target_cmap = mesh_->cell_map(true);

  cell_importer_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));

  const Epetra_Map& source_fmap = mesh_->face_map(false);
  const Epetra_Map& target_fmap = mesh_->face_map(true);

  face_importer_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
#endif

  // time control
  set_time(0.0, 1e-8);  // default parameters
  dT_desirable_ = dT;

  // miscalleneous
  mfd3d_method = FLOW_MFD3D_OPTIMIZED;  // will be changed (lipnikov@lanl.gov)
  verbosity = FLOW_VERBOSITY_HIGH;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Darcy_PK::~Darcy_PK()
{
  delete super_map_;
  delete solver;
  if (matrix == preconditioner) {
    delete matrix;
  } else {
    delete matrix;
    delete preconditioner;
  }
  delete bc_pressure;
  delete bc_head;
  delete bc_flux;
  delete bc_seepage;
}


/* ******************************************************************
* Extract information from Diffusion Problem parameter list.
****************************************************************** */
void Darcy_PK::InitPK(Matrix_MFD* matrix_, Matrix_MFD* preconditioner_)
{
  if (matrix_ == NULL)
    matrix = new Matrix_MFD(FS, *super_map_);
  else
    matrix = matrix_;

  if (preconditioner_ == NULL)
    preconditioner = matrix;
  else
    preconditioner = preconditioner_;

  // Create the solution vectors.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->CreateCellView(*solution));
  solution_faces = Teuchos::rcp(FS->CreateFaceView(*solution));

  solver = new AztecOO;
  solver->SetUserOperator(matrix);
  solver->SetPrecOperator(preconditioner);
  solver->SetAztecOption(AZ_solver, AZ_cg);

  // Get parameters from the flow parameter list.
  ProcessParameterList();

  // Process boundary data
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers.resize(nfaces, FLOW_BC_FACE_NULL);
  bc_values.resize(nfaces, 0.0);

  double time = FS->get_time();
  if (time >= 0.0) T_physics = time;

  time = T_physics;
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);
  ProcessBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *solution_faces, atm_pressure,
      bc_markers, bc_values);

  // Process other fundamental structures
  K.resize(ncells_owned);
  matrix->SetSymmetryProperty(true);
  matrix->SymbolicAssembleGlobalMatrices(*super_map_);

  // Allocate data for relative permeability (for consistency).
  Krel_cells = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(true)));
  Krel_faces = Teuchos::rcp(new Epetra_Vector(mesh_->face_map(true)));

  Krel_cells->PutScalar(1.0);
  Krel_faces->PutScalar(1.0);  // must go away (lipnikov@lanl.gov)

  // Preconditioner
  Teuchos::ParameterList ML_list = preconditioner_list_.sublist(preconditioner_name_sss_).sublist("ML Parameters");
  preconditioner->InitML_Preconditioner(ML_list);

  flow_status_ = FLOW_STATUS_INIT;
};


/* ******************************************************************
* Separate initialization of solver may be required for steady state
* and transient runs.       
****************************************************************** */
void Darcy_PK::InitSteadyState(double T0, double dT0)
{
  set_time(T0, dT0);
  dT_desirable_ = dT0;  // The desirable time step from now on.
  num_itrs_sss = 0;

  Epetra_Vector& pressure = FS->ref_pressure();
  *solution_cells = pressure;

  // initialize mass matrices
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < K.size(); c++) K[c] *= rho_ / mu_;
  matrix->CreateMFDmassMatrices(mfd3d_method, K);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
    double pokay = 100 * matrix->nokay() / double(ncells_owned);
    double ppassed = 100 * matrix->npassed() / double(ncells_owned);
    std::printf("Darcy PK: Successful plus passed matrices: %4.1f%% %4.1f%%\n", pokay, ppassed);   
  }

  flow_status_ = FLOW_STATUS_STEADY_STATE_INIT;
}


/* ******************************************************************
* Initialization analyzes status of matrix/preconditioner pair.      
****************************************************************** */
void Darcy_PK::InitTransient(double T0, double dT0)
{
  set_time(T0, dT0);
  dT_desirable_ = dT0;  // The desirable time step from now on.
  num_itrs_trs = 0;

  // initialize mass matrices
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < K.size(); c++) K[c] *= rho_ / mu_;
  matrix->CreateMFDmassMatrices(mfd3d_method, K);

  flow_status_ = FLOW_STATUS_TRANSIENT_STATE_INIT;
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
int Darcy_PK::AdvanceToSteadyState()
{
  solver->SetAztecOption(AZ_output, AZ_none);

  // calculate and assemble elemental stifness matrices
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  matrix->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix);
  matrix->ApplyBoundaryConditions(bc_markers, bc_values);
  matrix->AssembleGlobalMatrices();
  matrix->ComputeSchurComplement(bc_markers, bc_values);
  matrix->UpdateML_Preconditioner();

  rhs = matrix->rhs();
  Epetra_Vector b(*rhs);
  solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
  solver->SetLHS(&*solution);  // initial solution guess

  solver->Iterate(max_itrs_sss, convergence_tol_sss);
  num_itrs_sss = solver->NumIters();
  residual_sss = solver->TrueResidual();

  if (verbosity >= FLOW_VERBOSITY_HIGH && MyPID == 0) {
    std::cout << "Darcy solver performed " << num_itrs_sss << " iterations." << std::endl
              << "Norm of true residual = " << residual_sss << std::endl;
  }

  flow_status_ = FLOW_STATUS_STEADY_STATE_COMPLETE;
  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT. The boundary conditions are 
* calculated only once, during the initialization step.  
******************************************************************* */
int Darcy_PK::Advance(double dT_MPC) 
{
  dT = dT_MPC;
  double time = FS->get_time();
  if (time >= 0.0) T_physics = time;

  solver->SetAztecOption(AZ_output, AZ_none);

  // update boundary conditions and source terms
  time = T_physics;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);

  if (src_sink != NULL) src_sink->Compute(time);

  ProcessBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *solution_faces, atm_pressure,
      bc_markers, bc_values);

  // calculate and assemble elemental stifness matrices
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  matrix->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix);
  AddTimeDerivativeSpecificStorage(*solution_cells, dT, matrix);
  matrix->ApplyBoundaryConditions(bc_markers, bc_values);
  matrix->AssembleGlobalMatrices();
  matrix->ComputeSchurComplement(bc_markers, bc_values);
  matrix->UpdateML_Preconditioner();

  rhs = matrix->rhs();
  if (src_sink != NULL) AddSourceTerms(src_sink, *rhs);

  Epetra_Vector b(*rhs);
  solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
  solver->SetLHS(&*solution);  // initial solution guess

  solver->Iterate(max_itrs_sss, convergence_tol_sss);
  num_itrs_sss = solver->NumIters();
  residual_sss = solver->TrueResidual();

  num_itrs_trs++;

  flow_status_ = FLOW_STATUS_TRANSIENT_STATE_COMPLETE;
  return 0;
}


/* ******************************************************************
* Transfer data from the external flow state FS_MPC. MPC may request
* to populate the original state FS. 
****************************************************************** */
void Darcy_PK::CommitState(Teuchos::RCP<Flow_State> FS_MPC)
{
  FS_MPC->ref_pressure() = *solution_cells;

  // calculate darcy mass flux
  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  matrix->DeriveDarcyMassFlux(*solution, *face_importer_, flux);
  AddGravityFluxes_DarcyFlux(K, *Krel_cells, *Krel_faces, flux);
  for (int c = 0; c < nfaces_owned; c++) flux[c] /= rho_;

  // DEBUG
  // WriteGMVfile(FS_MPC);
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Darcy_PK::SetAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& vertical_permeability = FS->ref_vertical_permeability();
  const Epetra_Vector& horizontal_permeability = FS->ref_horizontal_permeability();

  for (int c = 0; c < K.size(); c++) {
    if (vertical_permeability[c] == horizontal_permeability[c]) {
      K[c].init(dim, 1);
      K[c](0, 0) = vertical_permeability[c];
    } else {
      K[c].init(dim, 2);
      for (int i = 0; i < dim-1; i++) K[c](i, i) = horizontal_permeability[c];
      K[c](dim-1, dim-1) = vertical_permeability[c];
    }
  }
}


/* ******************************************************************
* Adds time derivative to cell-based part of MFD algebraic system.
* Specific storage at the moment is 1.                                              
****************************************************************** */
void Darcy_PK::AddTimeDerivativeSpecificStorage(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix)
{
  const Epetra_Vector& specific_storage = FS->ref_specific_storage();

  std::vector<double>& Acc_cells = matrix->Acc_cells();
  std::vector<double>& Fc_cells = matrix->Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume / dT_prec * specific_storage[c];
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* A wrapper for a similar matrix call.  
****************************************************************** */
void Darcy_PK::DeriveDarcyVelocity(const Epetra_Vector& flux, Epetra_MultiVector& velocity)
{
  matrix->DeriveDarcyVelocity(flux, *face_importer_, velocity);
}


/* ******************************************************************
* Printing information about Flow status.                                                     
****************************************************************** */
void Darcy_PK::PrintStatistics() const
{
  if (!MyPID && verbosity > 0) {
    cout << "Flow PK:" << endl;
    cout << "    Verbosity level = " << verbosity << endl;
    cout << "    Enable internal tests = " << (internal_tests ? "yes" : "no")  << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

