/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Usage: Richards_PK FPK(parameter_list, state);
         FPK->InitPK();
         FPK->Initialize();  // optional
         FPK->InitSteadyState(T, dT);
         FPK->InitTransient(T, dT);
*/

#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "mfd3d_diffusion.hh"

#include "Matrix.hh"
#include "MatrixFactory.hh"
#include "Flow_BC_Factory.hh"
#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Simplest possible constructor: extracts lists and requires fields.
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S)
{
  S_ = S;

  mesh_ = S_->GetMesh();
  dim = mesh_->space_dimension();

  // We need the flow list
  Teuchos::ParameterList flow_list;
  if (glist.isSublist("Flow")) {
    flow_list = glist.sublist("Flow");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Richards Problem")) {
    rp_list_ = flow_list.sublist("Richards Problem");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Richards Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  // We also need iscaleneous sublists
  if (glist.isSublist("Preconditioners")) {
    preconditioner_list_ = glist.sublist("Preconditioners");
  } else {
    Errors::Message msg("Flow PK: input XML does not have <Preconditioners> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (glist.isSublist("Solvers")) {
    linear_operator_list_ = glist.sublist("Solvers");
  } else {
    Errors::Message msg("Flow PK: input XML does not have <Solvers> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  // for creating fields
  std::vector<std::string> names(2);
  names[0] = "cell"; 
  names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL; 
  locations[1] = AmanziMesh::FACE;

  std::vector<int> ndofs(2, 1);

  // require state variables for the Richards PK
  if (!S_->HasField("fluid_density")) {
    S_->RequireScalar("fluid_density", passwd_);
  }
  if (!S_->HasField("fluid_voscosity")) {
    S_->RequireScalar("fluid_viscosity", passwd_);
  }
  if (!S_->HasField("gravity")) {
    S_->RequireConstantVector("gravity", passwd_, dim);
  }

  if (!S_->HasField("pressure")) {
    S_->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }
  if (!S_->HasField("hydraulic_head")) {
    S_->RequireField("hydraulic_head", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);
  }

  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("water_saturation")) {
    S_->RequireField("water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("prev_water_saturation")) {
    S_->RequireField("prev_water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  if (!S_->HasField("darcy_velocity")) {
    S_->RequireField("darcy_velocity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, mesh_->space_dimension());
  }
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  delete bc_pressure;
  delete bc_flux;
  delete bc_head;
  delete bc_seepage;

  if (ti_specs != NULL) OutputTimeHistory(rp_list_, ti_specs->dT_history);
}


/* ******************************************************************
* Extract information from Richards Problem parameter list.
****************************************************************** */
void Richards_PK::InitPK()
{
  // Initialize miscalleneous default parameters.
  ti_specs = NULL;
  block_picard = 1;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  mfd3d_method_ = FLOW_MFD3D_OPTIMIZED;
  mfd3d_method_preconditioner_ = FLOW_MFD3D_OPTIMIZED;

  src_sink = NULL;
  src_sink_distribution = 0;

  // Initilize various common data depending on mesh and state.
  Flow_PK::Init();

  // Time control specific to this PK.
  ResetPKtimes(0.0, FLOW_INITIAL_DT);
  dT_desirable_ = dT;

  // Allocate memory for boundary data.
  bc_tuple zero = {0.0, 0.0};
  bc_values.resize(nfaces_wghost, zero);
  bc_model.resize(nfaces_wghost, 0);
  bc_submodel.resize(nfaces_wghost, 0);

  rainfall_factor.resize(nfaces_wghost, 1.0);

  // Process Native XML.
  ProcessParameterList(rp_list_);

  // Create water retention models.
  Teuchos::ParameterList& wrm_list = rp_list_.sublist("Water retention models");
  rel_perm = Teuchos::rcp(new RelativePermeability(mesh_, wrm_list));
  rel_perm->Init(atm_pressure_, S_);
  rel_perm->ProcessParameterList();
  rel_perm->PopulateMapC2MB();

  std::string krel_method_name = rp_list_.get<string>("relative permeability");
  rel_perm->ProcessStringRelativePermeability(krel_method_name);

  // Initialize times.
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // Initialize actions on boundary condtions. 
  ProcessShiftWaterTableList(rp_list_);

  time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(time);
  else
    bc_head->ComputeShift(time, shift_water_table_->Values());

  // Process other fundamental structures.
  K.resize(ncells_wghost);
  SetAbsolutePermeabilityTensor();
  rel_perm->CalculateKVectorUnit(K, gravity_);

  // Allocate memory for wells.
  const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false);
  if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  }

  // Select a proper matrix class. 
  Teuchos::ParameterList mlist;
  if (mfd3d_method_ == FLOW_FV_TPFA) {
    mlist.set<std::string>("matrix", "tpfa");
  } else {
    mlist.set<std::string>("matrix", "mfd");
  }

  MatrixFactory factory;
  matrix_ = factory.Create(S_, &K, rel_perm, mlist);
  preconditioner_ = factory.Create(S_, &K, rel_perm, mlist);

  matrix_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_MATRIX);
  preconditioner_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER);

  is_matrix_symmetric = SetSymmetryProperty();
  matrix_->SetSymmetryProperty(is_matrix_symmetric);
  matrix_->SymbolicAssemble();

  // Create the solution (pressure) vector and auxiliary vector for time history.
  solution = Teuchos::rcp(new CompositeVector(matrix_->DomainMap()));
  solution->PutScalar(0.0);
  
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // Initialize boundary and source data. 
  CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
  UpdateSourceBoundaryData(time, pressure);

  // Other quantatities: injected water mass
  mass_bc = 0.0;
}


/* ******************************************************************
* Initialization of auxiliary variables (lambda and two saturations).
* WARNING: Flow_PK may use complex initialization of the remaining 
* state variables.
****************************************************************** */
void Richards_PK::InitializeAuxiliaryData()
{
  // pressures
  CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
  const Epetra_MultiVector& p_cells = *pressure.ViewComponent("cell");
  Epetra_MultiVector& p_faces = *pressure.ViewComponent("face");

  DeriveFaceValuesFromCellValues(p_cells, p_faces);

  double time = T_physics;
  UpdateSourceBoundaryData(time, pressure);

  // saturations
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p_cells, ws);
  ws_prev = ws;
}


/* ******************************************************************
* Initial pressure is set to the pressure for fully saturated rock.
****************************************************************** */
void Richards_PK::InitializeSteadySaturated()
{ 
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "initializing with a saturated steady state..." << endl;
  }
  double T = S_->time();
  SolveFullySaturatedProblem(T, *solution, ti_specs->ls_specs_ini);
}


/* ******************************************************************
* Specific initialization of the initial pressure guess calculation.
****************************************************************** */
void Richards_PK::InitPicard(double T0)
{
  ti_specs = &ti_specs_igs_;
  if (ti_specs->error_control_options == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;
  }

  InitNextTI(T0, 0.0, ti_specs_igs_);

  // calculate initial guess: cleaning is required (lipnikov@lanl.gov)
  T_physics = ti_specs_igs_.T0;
  dT = ti_specs_igs_.dT0;
  AdvanceToSteadyState_Picard(ti_specs_igs_);

  const Epetra_MultiVector& p_cells = *S_->GetFieldData("pressure")->ViewComponent("cell");
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p_cells, ws);
  ws_prev = ws;
}


/* ******************************************************************
* Specific initialization of a steady state time integration phase.
****************************************************************** */
void Richards_PK::InitSteadyState(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(rp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_sss_;

  error_control_ = ti_specs->error_control_options;
  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
    ti_specs->error_control_options = error_control_;
  }

  InitNextTI(T0, dT0, ti_specs_sss_);
}


/* ******************************************************************
* Specific initialization of a transient time integration phase.  
* WARNING: Initialization of lambda is done in MPC and may be 
* erroneous in pure transient mode.
****************************************************************** */
void Richards_PK::InitTransient(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(rp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_trs_;

  error_control_ = ti_specs->error_control_options;
  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4
    ti_specs->error_control_options = error_control_;
  }

  InitNextTI(T0, dT0, ti_specs_trs_);
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Richards_PK::InitNextTI(double T0, double dT0, TI_Specs& ti_specs)
{
  ResetPKtimes(T0, dT0);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    LinearSolver_Specs& ls = ti_specs.ls_specs;
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << endl 
        << "****************************************" << endl
        << vo_->color("green") << "New TI phase: " << ti_specs.ti_method_name.c_str() << vo_->reset() << endl
        << "****************************************" << endl
        << " start T=" << T0 / FLOW_YEAR << " [y], dT=" << dT0 << " [sec]" << endl
        << " error control id=" << error_control_ << endl
        << " preconditioner for nonlinear solver: " << ti_specs.preconditioner_name.c_str() << endl
        << " sources distribution id=" << src_sink_distribution << endl;

    if (ti_specs.initialize_with_darcy) {
      LinearSolver_Specs& ls_ini = ti_specs.ls_specs_ini;
      *vo_->os() << " initial pressure solver: " << ls_ini.solver_name << endl;
      if (ti_specs.clip_saturation > 0.0) {
        *vo_->os() << "  clipping saturation at " << ti_specs.clip_saturation << " [-]" << endl;
      } else if (ti_specs.clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
        *vo_->os() << "  clipping pressure at " << ti_specs.clip_pressure << " [Pa]" << endl;
      }
    }
  }

  // set up new preconditioner
  Teuchos::ParameterList& tmp_list = preconditioner_list_.sublist(ti_specs.preconditioner_name);
  string mfd3d_method_name = tmp_list.get<string>("discretization method", "monotone mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_preconditioner_); 

  // preconditioner_->DestroyPreconditioner();
  preconditioner_->SetSymmetryProperty(is_matrix_symmetric);
  preconditioner_->SymbolicAssemble();
  preconditioner_->InitPreconditioner(ti_specs.preconditioner_name, preconditioner_list_);

  // set up new time integration or solver
  std::string ti_method_name(ti_specs.ti_method_name);

  if (ti_specs.ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList bdf1_list = rp_list_.sublist(ti_method_name).sublist("BDF1");
    if (! bdf1_list.isSublist("VerboseObject"))
        bdf1_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    bdf1_dae = Teuchos::rcp(new BDF1_TI<CompositeVector, CompositeVectorSpace>(*this, bdf1_list, solution));
  }

  // initialize mass matrices
  SetAbsolutePermeabilityTensor();
  for (int c = 0; c < ncells_wghost; c++) K[c] *= rho_ / mu_;

  matrix_->CreateMassMatrices(mfd3d_method_);
  preconditioner_->CreateMassMatrices(mfd3d_method_preconditioner_);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    int missed_tmp = missed_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
#endif

    int nokay = matrix_->nokay();
    int npassed = matrix_->npassed();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << " discretization method (prec): " << mfd3d_method_name.c_str() << endl
               << "  good and repaired matrices: " << nokay << " " << npassed << endl
               << " assigned default (no-flow) BC to " << missed_bc_faces_ << " faces" << endl << endl;
  }

  // Well modeling
  if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    CalculatePermeabilityFactorInWell();
  }

  // Initialize pressure (p and lambda components of solution and State).
  Epetra_MultiVector& pstate = *S_->GetFieldData("pressure", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& p = *solution->ViewComponent("cell");
  p = pstate;

  bool flag_face(false);
  if (solution->HasComponent("face")) {
    flag_face = true;
    *solution->ViewComponent("face") = *S_->GetFieldData("pressure")->ViewComponent("face");
  }

  if (ti_specs.initialize_with_darcy) {
    // Get a hydrostatic solution consistent with b.c.
    SolveFullySaturatedProblem(T0, *solution, ti_specs.ls_specs_ini);

    bool clip(false);
    if (ti_specs.clip_saturation > 0.0) {
      double pmin = FLOW_PRESSURE_ATMOSPHERIC;
      ClipHydrostaticPressure(pmin, ti_specs.clip_saturation, p);
      clip = true;
    } else if (ti_specs.clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
      ClipHydrostaticPressure(ti_specs.clip_pressure, p);
      clip = true;
    }
    pstate = p;

    if (flag_face && clip) {
      Epetra_MultiVector& lambda = *solution->ViewComponent("face", true);
      DeriveFaceValuesFromCellValues(p, lambda);
    }
  }


  // initialize saturation
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  DeriveSaturationFromPressure(pstate, ws);
 
  // re-initialize lambda (experimental)
  double Tp = T0 + dT0;
 if (flag_face && ti_specs.pressure_lambda_constraints) {
    EnforceConstraints(Tp, *solution);
  } else {
    CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
    UpdateSourceBoundaryData(Tp, *solution);
    rel_perm->Compute(pressure, bc_model, bc_values);
  }

  // nonlinear solver control options
  ti_specs.num_itrs = 0;
  block_picard = 0;

  CompositeVector& darcy_flux = *S_->GetFieldData("darcy_flux", passwd_);
  matrix_->CreateStiffnessMatricesRichards();
  matrix_->DeriveMassFlux(*solution, darcy_flux, bc_model, bc_values);

  Epetra_MultiVector& flux = *darcy_flux.ViewComponent("face", true);
  AddGravityFluxes_DarcyFlux(flux, *rel_perm);

  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;
}


/* ******************************************************************
* This routine avoid limitations of MPC by bumping up the time step.                                          
****************************************************************** */
double Richards_PK::CalculateFlowDt()
{
  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD && block_picard == 1) dT *= 1e+4;
  return dT;
}


/* ******************************************************************* 
* Performs one time step of size dT_MPC either for steady-state or 
* transient calculations.
* Warning: BDF2 and BDF1 will merge eventually.
******************************************************************* */
int Richards_PK::Advance(double dT_MPC)
{
  Teuchos::RCP<CompositeVector> cv = S_->GetFieldData("darcy_flux", passwd_);
  cv->ScatterMasterToGhosted("face");
  Epetra_MultiVector& flux = *cv->ViewComponent("face", true);

  dT = dT_MPC;
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // predict water mass change during time step
  time = T_physics;
  if (ti_specs->num_itrs == 0) {  // initialization
    // I do not know how to estimate du/dt, so it is set to zero.
    Teuchos::RCP<CompositeVector> udot = Teuchos::rcp(new CompositeVector(*solution));
    udot->PutScalar(0.0);

    if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
      bdf1_dae->SetInitialState(time, solution, udot);

    } else if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD) {
      AdvanceToSteadyState(time, dT_MPC);
      block_picard = 1;  // We will wait for transient initialization.
    }

    UpdatePreconditioner(time, solution, dT);
    ti_specs->num_itrs++;
  }

  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    while (bdf1_dae->TimeStep(dT, dTnext, solution)) {
      dT = dTnext;
    }
    bdf1_dae->CommitSolution(dT, solution);
    T_physics = bdf1_dae->time();
  }

  // Calculate time derivative and 2nd-order solution approximation.
  // Estimate of a time step multiplier overrides the above estimate.
  if (ti_specs->dT_method == FLOW_DT_ADAPTIVE) {
    const Epetra_MultiVector& pressure = *S_->GetFieldData("pressure")->ViewComponent("face", false);
    Epetra_MultiVector& p = *solution->ViewComponent("face", false);

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = (p[0][c] - pressure[0][c]) / dT; 
      p[0][c] = pressure[0][c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dT / 2;
    }

    double err, dTfactor;
    err = AdaptiveTimeStepEstimate(&dTfactor);
    if (err > 0.0) throw 1000;  // fix (lipnikov@lan.gov)
    dTnext = std::min(dT_MPC * dTfactor, ti_specs->dTmax);
  }

  dt_tuple times(time, dT_MPC);
  ti_specs->dT_history.push_back(times);

  ti_specs->num_itrs++;

  return 0;
}


/* ******************************************************************
* Transfer part of the internal data needed by transport to the 
* flow state FS_MPC. MPC may request to populate the original FS.
* The consistency condition is improved by adjusting saturation while
* preserving its LED property.
****************************************************************** */
void Richards_PK::CommitState(Teuchos::RCP<State> S)
{
  // copy solution to State
  CompositeVector& pressure = *S->GetFieldData("pressure", passwd_);
  *pressure.ViewComponent("cell") = *solution->ViewComponent("cell");

  if (solution->HasComponent("face")) {
    *pressure.ViewComponent("face") = *solution->ViewComponent("face");
  }

  // ws -> ws_prev
  Epetra_MultiVector& ws = *S->GetFieldData("water_saturation", passwd_)->ViewComponent("cell", false);
  Epetra_MultiVector& ws_prev = *S->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell", false);
  ws_prev = ws;

  // calculate new water saturation
  const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell", false);
  DeriveSaturationFromPressure(p, ws);

  // calculate Darcy flux as diffusive part + advective part.
  CompositeVector& darcy_flux = *S->GetFieldData("darcy_flux", passwd_);
  matrix_->CreateStiffnessMatricesRichards();
  matrix_->DeriveMassFlux(*solution, darcy_flux, bc_model, bc_values);

  Epetra_MultiVector& flux = *darcy_flux.ViewComponent("face", true);
  AddGravityFluxes_DarcyFlux(flux, *rel_perm);

  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  // update mass balance
  // ImproveAlgebraicConsistency(ws_prev, ws);
  
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Epetra_MultiVector& phi = *S_->GetFieldData("porosity", passwd_)->ViewComponent("cell", false);
    double mass_bc_dT = WaterVolumeChangePerSecond(bc_model, flux) * rho_ * dT;

    mass_amanzi = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      mass_amanzi += ws[0][c] * rho_ * phi[0][c] * mesh_->cell_volume(c);
    }

    double mass_amanzi_tmp = mass_amanzi, mass_bc_tmp = mass_bc_dT;
    mesh_->get_comm()->SumAll(&mass_amanzi_tmp, &mass_amanzi, 1);
    mesh_->get_comm()->SumAll(&mass_bc_tmp, &mass_bc_dT, 1);

    mass_bc += mass_bc_dT;

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "water mass=" << mass_amanzi 
                 << " [kg], total boundary flux=" << mass_bc << " [kg]" << endl;
  }

  dT = dTnext;
}


/* ******************************************************************
* Estimate du/dt from the pressure equation, a du/dt = g - A*u.
* There is no meaniful way to estimate du/dt due to the abrupt change 
* of function a. The routine is not used.
****************************************************************** */
/*
double Richards_PK::ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot)
{
  AssembleMatrixMFD(u, T);
  double norm_udot = matrix_->ComputeNegativeResidual(u, udot);

  Epetra_Vector* udot_cells = FS->CreateCellView(udot);
  const Epetra_Vector& phi = FS->ref_porosity();
  Epetra_Vector dSdP(mesh_->cell_map(false));
  rel_perm->DerivedSdP(u, dSdP);
 
  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    if (dSdP[c] > 0.0) (*udot_cells)[c] /= volume * dSdP[c] * phi[c] * rho_;
  }

  Epetra_Vector* udot_faces = FS->CreateFaceView(udot);
  DeriveFaceValuesFromCellValues(*udot_cells, *udot_faces);

  return norm_udot;
}
*/


/* ******************************************************************
* Identify type of the matrix.
****************************************************************** */
bool Richards_PK::SetSymmetryProperty()
{
  int method = rel_perm->method();
  bool sym = false;
  if ((method == FLOW_RELATIVE_PERM_CENTERED) ||
      (method == FLOW_RELATIVE_PERM_AMANZI)) sym = true;
  return sym;
}


/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations. 
* WARNING: it is implmented for transient phase only.
****************************************************************** */
double Richards_PK::AdaptiveTimeStepEstimate(double* dTfactor)
{
  Epetra_MultiVector& p = *solution->ViewComponent("cell");

  double tol, error, error_max = 0.0;
  double dTfactor_cell;

  *dTfactor = 100.0;
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs((*pdot_cells)[c] - (*pdot_cells_prev)[c]) * dT / 2;
    tol = ti_specs_trs_.rtol * fabs(p[0][c]) + ti_specs_trs_.atol;

    dTfactor_cell = sqrt(tol / std::max(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dTfactor = std::min(*dTfactor, dTfactor_cell);

    error_max = std::max(error_max, error - tol);
  }

  *dTfactor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dTfactor = std::min(*dTfactor, FLOW_DT_ADAPTIVE_INCREASE);
  *dTfactor = std::max(*dTfactor, FLOW_DT_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
    double dT_tmp = *dTfactor;
    solution->Comm().MinAll(&dT_tmp, dTfactor, 1);  // find the global minimum
 
    double error_tmp = error_max;
    solution->Comm().MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif
  return error_max;
}


/* ******************************************************************
* Saturation should be in exact balance with Darcy fluxes in order to
* have local extrema dimishing (LED) property for concentrations. 
* WARNING: we can enforce it strictly only in some cells.
****************************************************************** */
void Richards_PK::ImproveAlgebraicConsistency(const Epetra_Vector& ws_prev, Epetra_Vector& ws)
{
  // create a disctributed flux and water_saturation vectors
  S_->GetFieldData("darcy_flux", passwd_)->ScatterMasterToGhosted("face");
  S_->GetFieldData("wateri_saturation", passwd_)->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);
  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell", false);

  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    // calculate min/max values
    double wsmin, wsmax;
    wsmin = wsmax = ws[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int c2 = mfd.cell_get_face_adj_cell(c, f);
      if (c2 >= 0) {
        wsmin = std::min(wsmin, ws[c2]);
        wsmax = std::max(wsmax, ws[c2]);
      }
    }

    // predict new saturation
    ws[c] = ws_prev[c];
    double factor = dT / (phi[0][c] * mesh_->cell_volume(c));
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      ws[c] -= factor * flux[0][f] * dirs[n]; 
    }

    // limit new saturation
    ws[c] = std::max(ws[c], wsmin);
    ws[c] = std::min(ws[c], wsmax);
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

