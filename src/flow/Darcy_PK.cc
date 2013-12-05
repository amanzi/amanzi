/*
This is the flow component of the Amanzi code.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_Import.h"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "exceptions.hh"

#include "mfd3d_diffusion.hh"
#include "tensor.hh"
#include "LinearOperatorFactory.hh"

#include "Darcy_PK.hh"
#include "FlowDefs.hh"
#include "Flow_SourceFactory.hh"
#include "Matrix_MFD.hh"
#include "Matrix_TPFA.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* each variable initialization
****************************************************************** */
Darcy_PK::Darcy_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S)
{
  // initialize pointers (Do we need smart pointers here? lipnikov@lanl.gov)
  bc_pressure = NULL; 
  bc_head = NULL;
  bc_flux = NULL;
  bc_seepage = NULL; 
  src_sink = NULL;

  Flow_PK::Init(glist, S);  // sets up default parameters

  // extract important sublists
  Teuchos::ParameterList flow_list;
  if (glist.isSublist("Flow")) {
    flow_list = glist.sublist("Flow");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Darcy Problem")) {
    dp_list_ = flow_list.sublist("Darcy Problem");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Darcy Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  // time control
  ResetPKtimes(0.0, FLOW_INITIAL_DT);
  dT_desirable_ = dT;

  // miscalleneous
  ti_specs = NULL;
  mfd3d_method_ = FLOW_MFD3D_OPTIMIZED;
  src_sink = NULL;
  src_sink_distribution = 0;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Darcy_PK::~Darcy_PK()
{
  delete bc_pressure;
  delete bc_head;
  delete bc_flux;
  delete bc_seepage;

  if (ti_specs != NULL) OutputTimeHistory(dp_list_, ti_specs->dT_history);
}


/* ******************************************************************
* Extract information from Diffusion Problem parameter list.
****************************************************************** */
void Darcy_PK::InitPK()
{
  // Allocate memory for boundary data. It must go first.
  bc_tuple zero = {0.0, 0.0};
  bc_values.resize(nfaces_wghost, zero);
  bc_model.resize(nfaces_wghost, 0);
  bc_submodel.resize(nfaces_wghost, 0);

  rainfall_factor.resize(nfaces_owned, 1.0);

  // Read flow list and populate various structures. 
  ProcessParameterList(dp_list_);

  // Select a proper matrix class. No optionos at the moment.
  matrix_ = Teuchos::rcp(new Matrix_MFD(mesh_));
  matrix_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_MATRIX);
  matrix_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER);

  // Create auxiliary data for time history.
  solution = Teuchos::rcp(new CompositeVector(*(S_->GetFieldData("pressure"))));

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  
  // Initialize times.
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // Initialize actions on boundary condtions. 
  ProcessShiftWaterTableList(dp_list_);

  time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(time);
  else
    bc_head->ComputeShift(time, shift_water_table_->Values());

  const CompositeVector& pressure = *S_->GetFieldData("pressure");
  ComputeBCs(pressure);

  // Process other fundamental structures
  K.resize(ncells_owned);
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();

  // Allocate memory for wells
  if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  }
};


/* ******************************************************************
* Initialization of auxiliary variables (lambda and two saturations).
* WARNING: Flow_PK may use complex initialization of the remaining 
* state variables.
****************************************************************** */
void Darcy_PK::InitializeAuxiliaryData()
{
  // pressures (lambda is not important when solver is very accurate)
  CompositeVector& cv = *S_->GetFieldData("pressure", passwd_);
  const Epetra_MultiVector& pressure = *(cv.ViewComponent("cell"));
  Epetra_MultiVector& lambda = *(cv.ViewComponent("face"));

  DeriveFaceValuesFromCellValues(pressure, lambda);

  // saturations
  if (!S_->GetField("water_saturation", passwd_)->initialized()) {
    S_->GetFieldData("water_saturation", passwd_)->PutScalar(1.0);
    S_->GetField("water_saturation", passwd_)->set_initialized();
  }
  if (!S_->GetField("prev_water_saturation", passwd_)->initialized()) {
    S_->GetFieldData("prev_water_saturation", passwd_)->PutScalar(1.0);
    S_->GetField("prev_water_saturation", passwd_)->set_initialized();
  }

  // miscalleneous
  UpdateSpecificYield();
}


/* ******************************************************************
* Wrapper for a steady-state solver
****************************************************************** */
void Darcy_PK::InitializeSteadySaturated()
{ 
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "initializing with a saturated steady state..." << endl;
  }
  double T = S_->time();
  SolveFullySaturatedProblem(T, *solution);
}


/* ******************************************************************
* Specific initialization of a steady state time integration phase.
* WARNING: now it is equivalent to transient phase.
****************************************************************** */
void Darcy_PK::InitSteadyState(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(dp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_sss_;

  InitNextTI(T0, dT0, ti_specs_sss_);

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;  // usually 1e-4;
}


/* ******************************************************************
* Specific initialization of a transient time integration phase.  
****************************************************************** */
void Darcy_PK::InitTransient(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(dp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_trs_;

  InitNextTI(T0, dT0, ti_specs_trs_);

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;  // usually 1e-4

  // DEBUG
  // SolveFullySaturatedProblem(0.0, *solution);
  // CommitState(FS); WriteGMVfile(FS); exit(0);
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Darcy_PK::InitNextTI(double T0, double dT0, TI_Specs ti_specs)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    LinearSolver_Specs& ls_specs = ti_specs.ls_specs;

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "****************************************" << endl
                 << "TI phase: " << ti_specs.ti_method_name.c_str() << endl
                 << "****************************************" << endl
                 << "  start T=" << T0 / FLOW_YEAR << " [y], dT=" << dT0 << " [sec]" << endl
                 << "  time stepping id=" << ti_specs.dT_method << endl
                 << "  sources distribution id=" << src_sink_distribution << endl
                 << "  linear solver: ||r||<" << ls_specs.convergence_tol << " #itr<" << ls_specs.max_itrs << endl
                 << "  preconditioner: " << ti_specs.preconditioner_name.c_str() << endl;
    if (ti_specs.initialize_with_darcy) {
      *(vo_->os()) << "  initial pressure guess: \"saturated solution\"" << endl;
    }
  }

  // set up new preconditioner
  // matrix_->DestroyPreconditionerNew();
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->InitPreconditioner(ti_specs.preconditioner_name, preconditioner_list_);

  // set up initial guess for solution
  Epetra_MultiVector& p = *S_->GetFieldData("pressure", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& p_cells = *solution->ViewComponent("cell");
  p_cells = p;

  ResetPKtimes(T0, dT0);
  dT_desirable_ = dT0;  // The minimum desirable time step from now on.
  ti_specs.num_itrs = 0;

  // initialize mass matrices
  double rho = *S_->GetScalarData("fluid_density");
  double mu = *S_->GetScalarData("fluid_viscosity");

  SetAbsolutePermeabilityTensor();
  for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;

  matrix_->CreateMFDmassMatrices(mfd3d_method_, K);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    int nokay = matrix_->nokay();
    int npassed = matrix_->npassed();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "  good and repaired matrices: " << nokay << " " << npassed << endl;
  }

  // Well modeling (one-time call)
  if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    CalculatePermeabilityFactorInWell();
  }

  // Initialize source
  if (src_sink != NULL) {
    if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(T0, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(T0, NULL);
    }
  }
  
  // make initial guess consistent with boundary conditions
  if (ti_specs.initialize_with_darcy) {
    ti_specs.initialize_with_darcy = false;

    Epetra_MultiVector& p_faces = *solution->ViewComponent("face", true);
    DeriveFaceValuesFromCellValues(p, p_faces);

    SolveFullySaturatedProblem(T0, *solution);
    p = p_cells;
  }
}


/* ******************************************************************
* Wrapper for a steady-state solver
****************************************************************** */
int Darcy_PK::AdvanceToSteadyState(double T0, double dT0)
{ 
  ti_specs = &ti_specs_sss_;
  SolveFullySaturatedProblem(T0, *solution);
  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT. The boundary conditions are 
* calculated only once, during the initialization step.  
******************************************************************* */
int Darcy_PK::Advance(double dT_MPC) 
{
  dT = dT_MPC;
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // update boundary conditions and source terms
  time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(time);
  else
    bc_head->ComputeShift(time, shift_water_table_->Values());

  if (src_sink != NULL) {
    if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(time, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(time, NULL);
    }
  }

  ComputeBCs(*solution);

  // calculate and assemble elemental stifness matrices
  Epetra_MultiVector& p_cells = *solution->ViewComponent("cell");

  matrix_->CreateMFDstiffnessMatrices();
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(&*matrix_);
  AddTimeDerivativeSpecificStorage(p_cells, dT, &*matrix_);
  AddTimeDerivativeSpecificYield(p_cells, dT, &*matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->AssembleSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();

  CompositeVector& rhs = *matrix_->rhs();
  if (src_sink != NULL) AddSourceTerms(rhs);

  // create linear solver
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AmanziSolvers::LinearOperatorFactory<Matrix_MFD, CompositeVector, Epetra_BlockMap> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, CompositeVector, Epetra_BlockMap> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);
  solver->ApplyInverse(rhs, *solution);

  ti_specs->num_itrs++;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "pressure solver (" << solver->name()
                 << "): ||r||=" << residual << " itr=" << num_itrs << endl;
  }

  // calculate time derivative and 2nd-order solution approximation
  if (ti_specs->dT_method == FLOW_DT_ADAPTIVE) {
    const Epetra_MultiVector& p = *S_->GetFieldData("pressure")->ViewComponent("cell");  // pressure at t^n
    Epetra_MultiVector& p_cells = *solution->ViewComponent("cell");  // pressure at t^{n+1}

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = (p_cells[0][c] - p[0][c]) / dT; 
      p_cells[0][c] = p[0][c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dT / 2;
    }
  }

  // estimate time multiplier
  if (ti_specs->dT_method == FLOW_DT_ADAPTIVE) {
    double err, dTfactor;
    err = ErrorEstimate(&dTfactor);
    if (err > 0.0) throw 1000;  // fix (lipnikov@lan.gov)
    dT_desirable_ = std::min(dT_MPC * dTfactor, ti_specs->dTmax);
  } else {
    dT_desirable_ = std::min(dT_desirable_ * ti_specs->dTfactor, ti_specs->dTmax);
  }

  dt_tuple times(time, dT_MPC);
  ti_specs->dT_history.push_back(times);

  return 0;
}


/* ******************************************************************
* Transfer data from the external flow state FS_MPC. MPC may request
* to populate the original state FS. 
****************************************************************** */
void Darcy_PK::CommitState(Teuchos::RCP<State> S)
{
  CompositeVector& p = *S_->GetFieldData("pressure", passwd_);
  p = *solution;

  // calculate darcy mass flux
  double rho = *S_->GetScalarData("fluid_density");
  Epetra_MultiVector& flux = *S_->GetFieldData("pressure", passwd_)->ViewComponent("face", true);

  matrix_->CreateMFDstiffnessMatrices();
  matrix_->DeriveDarcyMassFlux(*solution, bc_model, bc_values, flux);
  AddGravityFluxes_DarcyFlux(flux);

  for (int c = 0; c < nfaces_owned; c++) flux[0][c] /= rho;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  // DEBUG
  // WriteGMVfile(FS_MPC);
}


/* ******************************************************************
* Adds time derivative related to specific stroage to cell-based 
* part of MFD algebraic system. 
****************************************************************** */
void Darcy_PK::AddTimeDerivativeSpecificStorage(
    Epetra_MultiVector& p_cells, double dTp, Matrix_MFD* matrix_operator)
{
  double g = fabs(gravity_[dim - 1]);
  const Epetra_MultiVector& 
      specific_storage = *S_->GetFieldData("specific_storage")->ViewComponent("cell");

  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * specific_storage[0][c] / (g * dTp);
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * p_cells[0][c];
  }
}


/* ******************************************************************
* Add area/length factor to specific yield. 
****************************************************************** */
void Darcy_PK::UpdateSpecificYield()
{
  // populate ghost cells
  S_->GetFieldData("specific_yield", passwd_)->ScatterMasterToGhosted();
  const Epetra_MultiVector& 
      specific_yield = *S_->GetFieldData("specific_yield", passwd_)->ViewComponent("cell", true);

  WhetStone::MFD3D_Diffusion mfd3d(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int negative_yield = 0;
  for (int c = 0; c < ncells_owned; c++) {
    if (specific_yield[0][c] > 0.0) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      double area = 0.0;
      int nfaces = faces.size();
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        int c2 = mfd3d.cell_get_face_adj_cell(c, f);

        if (c2 >= 0) {
          if (specific_yield[0][c2] <= 0.0)  // cell in the fully saturated layer
            area -= (mesh_->face_normal(f))[dim - 1] * dirs[n];
        }
      }
      specific_yield[0][c] *= area;
      if (area <= 0.0) negative_yield++;
    }
  }

#ifdef HAVE_MPI
  int negative_yield_tmp = negative_yield;
  mesh_->get_comm()->MaxAll(&negative_yield_tmp, &negative_yield, 1);
#endif
  if (negative_yield > 0) {
    Errors::Message msg;
    msg << "Flow PK: configuration of the yield region leads to negative yield interfaces.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Adds time derivative related to specific yiled to cell-based part 
* of MFD algebraic system. Area factor is alreafy inside Sy. 
****************************************************************** */
void Darcy_PK::AddTimeDerivativeSpecificYield(
    Epetra_MultiVector& p_cells, double dTp, Matrix_MFD* matrix_operator)
{
  double g = fabs(gravity_[dim - 1]);
  const Epetra_MultiVector& 
      specific_yield = *S_->GetFieldData("specific_yield")->ViewComponent("cell");

  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double factor = specific_yield[0][c] / (g * dTp);
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * p_cells[0][c];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

