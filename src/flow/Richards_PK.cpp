/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

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
  ti_bdf2 = new Interface_BDF2(this, rme_list);
  itrs_bdf2 = 0;

  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(rp_list->sublist("Time integrator")));
  bdf2_dae = new BDF2::Dae(*ti_bdf2, *super_map_);
  bdf2_dae->setParameterList(bdf2_list);

  // Allocate data for relative permeability
  Krel_cells = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(true)));
  if (upwind_Krel) {
    const Epetra_Map& fmap = mesh_->face_map(true);

    Krel_faces = Teuchos::rcp(new Epetra_Vector(fmap));
    upwind_cell = Teuchos::rcp(new Epetra_IntVector(fmap));
    downwind_cell = Teuchos::rcp(new Epetra_IntVector(fmap));

    identify_upwind_cells(*upwind_cell, *downwind_cell);
  }
}


/* ******************************************************************
*  Wrapper for advance to steady-state routines.                                                    
****************************************************************** */
int Richards_PK::advance_to_steady_state()
{
  if (method_sss == FLOW_STEADY_STATE_PICARD) {
    return advanceSteadyState_Picard();
  } else if (method_sss == FLOW_STEADY_STATE_BACKWARD_EULER) {
    return advanceSteadyState_BackwardEuler();
  }
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advanceSteadyState_Picard()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;
  Epetra_Vector  residual(*solution);

  solver->SetAztecOption(AZ_output, AZ_none);
  int itrs = 0;
  double L2norm, L2error = 1.0, L2error_prev = 1.0;
  double relaxation = 0.0;

  while (L2error > err_tol_sss && itrs < max_itrs_sss) {
    // work-around limited support for tensors
    setAbsolutePermeabilityTensor(K);
    if (upwind_Krel) {
      identify_upwind_cells(*upwind_cell, *downwind_cell);
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
    } else { 
      calculateRelativePermeability(*solution_cells);
    }
    for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;

    matrix->createMFDstiffnessMatrices(K);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    // check for convergence of non-linear residual
    rhs = matrix->get_rhs();
    L2error = matrix->computeResidual(solution_new, residual);
    (*rhs).Norm2(&L2norm);
    L2error /= L2norm;

    if (L2error > L2error_prev) relaxation = std::min(0.95, relaxation * 1.5); 
    else relaxation = std::max(0.05, relaxation * 0.9);  

    // call AztecOO solver
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, err_tol);
    num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();
cout << itrs << " res=" << L2error << "  cg# " << num_itrs << "  relax=" << relaxation << endl;

    for (int c=0; c<number_owned_cells; c++) {
      solution_new[c] = relaxation * solution_old[c] + (1.0 - relaxation) * solution_new[c];
      solution_old[c] = solution_new[c];
    }

    L2error_prev = L2error;
    itrs++;
  }    
for (int c=0; c<K.size(); c+=4) cout << solution_new[c] << endl; 
  
  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->createMFDstiffnessMatrices(K);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyFlux(*solution, *face_importer_, darcy_flux);
  addGravityFluxes(darcy_flux);

  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advanceSteadyState_BackwardEuler()
{
  solver->SetAztecOption(AZ_output, AZ_none);

  int itrs = 0;
  double T = 0.0, Tend = 1000000.0;
  dT = 1000.0;

  while (T < Tend) {
    // work-around limited support for tensors
    setAbsolutePermeabilityTensor(K);
    calculateRelativePermeability(*solution_cells);
    for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;

    matrix->createMFDstiffnessMatrices(K);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(matrix);
    addTimeDerivative_MFD(*solution_cells, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    // call AztecOO solver
    rhs = matrix->get_rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, err_tol);
    num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();
cout << itrs << " res=" << residual << "  cg=" << num_itrs << endl;

    T += dT;
    itrs++;
  }
for (int c=0; c<K.size(); c+=4) cout << (*solution_cells)[c] << endl; 

  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT. 
******************************************************************* */
int Richards_PK::advance(double dT) 
{
  solver->SetAztecOption(AZ_output, AZ_none);

  // work-around limited support for tensors
  setAbsolutePermeabilityTensor(K);
  calculateRelativePermeability(*solution_cells);

  for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;

  matrix->createMFDstiffnessMatrices(K);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(matrix);
  addTimeDerivative_MFD(*solution_cells, matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();
  matrix->computeSchurComplement(bc_markers, bc_values);
  matrix->update_ML_preconditioner();

  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;
  double dTnext;

  if (itrs_bdf2 == 0) {  // initialization
    Epetra_Vector pdot(*super_map_);
    computePDot(time, *solution, pdot);
    bdf2_dae->set_initial_state(time, *solution, pdot);

    int errc;
    ti_bdf2->update_precon(time, *solution, dT, errc);
  }

  bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
  bdf2_dae->commit_solution(dT, *solution);
  bdf2_dae->write_bdf2_stepping_statistics();

  itrs_bdf2++;
  return 0;
}


/* ******************************************************************
* Estimate dp/dt from the pressure equations.                                               
****************************************************************** */
void Richards_PK::computePDot(const double T, const Epetra_Vector& p, Epetra_Vector &pdot)
{
  double norm_pdot;
  norm_pdot = matrix->computeResidual(p, pdot);

  Epetra_Vector *pdot_face = FS->create_face_view(pdot);
  pdot_face->PutScalar(0.0);
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

    calculateGravityFluxes(c, K, *upwind_cell, gravity_flux, upwind_Krel);
    
    Epetra_SerialDenseVector& Ff = matrix->get_Ff_cells()[c];
    for (int n=0; n<nfaces; n++) Ff[n] += gravity_flux[n] * dirs[n];
  }
}


/* ******************************************************************
* Adds .                                               
****************************************************************** */
void Richards_PK::addTimeDerivative_MFD(Epetra_Vector& pressure_cells, Matrix_MFD* matrix)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  derivedSdP(pressure_cells, dSdP);
 
  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  std::vector<double>& Fc_cells = matrix->get_Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = dSdP[c] * phi[c] * volume / dT;
if (factor < 0) cout << "HERE:" << dSdP[c] << endl;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

