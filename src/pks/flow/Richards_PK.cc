/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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

#include "boost/math/tools/roots.hpp"
#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "exceptions.hh"
#include "Mesh.hh"
#include "mfd3d_diffusion.hh"
#include "OperatorDiffusionFactory.hh"
#include "Point.hh"
#include "UpwindFactory.hh"
#include "XMLParameterListWriter.hh"

#include "darcy_velocity_evaluator.hh"
#include "Flow_BC_Factory.hh"
#include "primary_variable_field_evaluator.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Simplest possible constructor: extracts lists and requires fields.
****************************************************************** */
Richards_PK::Richards_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const std::string& pk_list_name,
                         Teuchos::RCP<State> S) :
    glist_(glist),
    Flow_PK()
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(pk_list, pk_list_name, true);
  rp_list_ = Teuchos::sublist(flow_list, "Richards problem", true);
  
  new_mpc_driver = glist->get<bool>("new mpc driver", false);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);

  if (rp_list_->isSublist("time integrator")) {
    ti_list_ = rp_list_->sublist("time integrator");
  } 
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void Richards_PK::Setup()
{
  mesh_ = S_->GetMesh();
  dim = mesh_->space_dimension();

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
  if (!S_->HasField("fluid_viscosity")) {
    S_->RequireScalar("fluid_viscosity", passwd_);
  }
  if (!S_->HasField("gravity")) {
    S_->RequireConstantVector("gravity", passwd_, dim);  // state resets ownerships.
  }

  if (!S_->HasField("pressure")) {
    S_->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "pressure");
    pressure_eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("pressure", pressure_eval);
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

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "darcy_flux");
    darcy_flux_eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("darcy_flux", darcy_flux_eval);
  }

  if (!S_->HasField("darcy_velocity")) {
    S_->RequireField("darcy_velocity", "darcy_velocity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist;
    Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S_->SetFieldEvaluator("darcy_velocity", eval);
  }
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  if (bc_pressure != NULL) delete bc_pressure;
  if (bc_flux != NULL) delete bc_flux;
  if (bc_head != NULL) delete bc_head;
  if (bc_seepage != NULL) delete bc_seepage;

  if (ti_specs != NULL) {
    OutputTimeHistory(*rp_list_, ti_specs->dT_history);
  }

  if (src_sink != NULL) delete src_sink;
  if (vo_ != NULL) delete vo_;

#ifdef ENABLE_NATIVE_XML_OUTPUT
  /*
  if (glist_->sublist("Analysis").get<bool>("print unused parameters", false)) {
    std::cout << "\n***** Unused XML parameters *****" << std::endl;
    Teuchos::Amanzi_XMLParameterListWriter out; 
    out.unused(*rp_list_, std::cout);
    std::cout << std::endl;
  }
  */
#endif
}


/* ******************************************************************
* Extract information from Richards Problem parameter list. It is 
* broken into a few pieces that can be reused in InitXXX routines.
* Cleaning will be done after switch to 
****************************************************************** */
void Richards_PK::Initialize()
{
  // Initialize miscalleneous default parameters.
  ti_specs = NULL;
  block_picard = 1;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  src_sink = NULL;
  src_sink_distribution = 0;

  // Initilize various common data depending on mesh and state.
  Flow_PK::Initialize();

  // Time control specific to this PK.
  ResetPKtimes(0.0, FLOW_INITIAL_DT);
  dT_desirable_ = dT;

  // Allocate memory for boundary data.
  bc_model.resize(nfaces_wghost, 0);
  bc_submodel.resize(nfaces_wghost, 0);
  bc_value.resize(nfaces_wghost, 0.0);
  bc_mixed.resize(nfaces_wghost, 0.0);
  op_bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  rainfall_factor.resize(nfaces_wghost, 1.0);

  // create verbosity object
  Teuchos::ParameterList vlist;
  vlist.sublist("VerboseObject") = rp_list_->sublist("VerboseObject");
  vo_ = new VerboseObject("FlowPK::Richards", vlist); 

  // Process Native XML.
  ProcessParameterList(*rp_list_);

  // Create water retention models.
  Teuchos::ParameterList& wrm_list = rp_list_->sublist("water retention models");
  rel_perm_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_->Init(atm_pressure_, wrm_list);

  std::string krel_method_name = rp_list_->get<std::string>("relative permeability");
  rel_perm_->ProcessStringRelativePermeability(krel_method_name);

  // parameter which defines when update direction of update
  Teuchos::ParameterList upw_list = rp_list_->sublist("operators")
                                             .sublist("diffusion operator")
                                             .sublist("upwind");
  Operators::UpwindFactory<RelativePermeability> upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, rel_perm_, upw_list);

  std::string upw_upd = rp_list_->get<std::string>("upwind update", "every timestep");
  if (upw_upd == "every nonlinear iteration") update_upwind = FLOW_UPWIND_UPDATE_ITERATION;
  else update_upwind = FLOW_UPWIND_UPDATE_TIMESTEP;  

  // Initialize times.
  double T1 = S_->time(), T0 = T1 - dT;
  if (T1 >= 0.0) T_physics = T1;

  // Initialize actions on boundary condtions. 
  ProcessShiftWaterTableList(*rp_list_);

  T1 = T_physics;
  bc_pressure->Compute(T1);
  bc_flux->Compute(T1);
  bc_seepage->Compute(T1);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(T1);
  else
    bc_head->ComputeShift(T1, shift_water_table_->Values());

  // Process other fundamental structures.
  K.resize(ncells_owned);
  SetAbsolutePermeabilityTensor();

  // Allocate memory for wells.
  if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(true)));
  }

  // Select a proper matrix class. 
  const Teuchos::ParameterList& tmp_list = rp_list_->sublist("operators")
                                                    .sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  std::string name = rp_list_->get<std::string>("relative permeability");
  int upw_id(Operators::OPERATOR_UPWIND_FLUX);
  std::string upw_method("none");
  if (name == "upwind: darcy velocity") {
    upw_method = "standard";
  } else if (name == "upwind: gravity") {
    upw_method = "standard";
    upw_id = Operators::OPERATOR_UPWIND_CONSTANT_VECTOR;
  } else if (name == "upwind: artificial diffusion") {
    upw_method = "amanzi: artificial diffusion";
    oplist_pc.set<std::string>("upwind method", "artificial diffusion");
  } else if (name == "upwind: amanzi") {
    upw_method = "divk";
  } else if (name == "other: arithmetic average") {
    upw_method = "standard";
    upw_id = Operators::OPERATOR_UPWIND_ARITHMETIC_AVERAGE;
  }
  oplist_matrix.set<std::string>("upwind method", upw_method);
  oplist_pc.set<std::string>("upwind method", upw_method);

  Operators::OperatorDiffusionFactory opfactory;
  op_matrix_diff_ = opfactory.Create(mesh_, op_bc_, oplist_matrix, gravity_, upw_id);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_preconditioner_diff_ = opfactory.Create(mesh_, op_bc_, oplist_pc, gravity_, upw_id);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_preconditioner_));

  // Create the solution (pressure) vector and auxiliary vector for time history.
  solution = Teuchos::rcp(new CompositeVector(op_matrix_->DomainMap()));
  
  const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // Initialize boundary and source data. 
  CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
  UpdateSourceBoundaryData(T0, T1, pressure);

  // Initialize two fields for upwind operators.
  InitializeUpwind_();

  // Other quantatities: injected water mass
  mass_bc = 0.0;
  seepage_mass_ = 0.0;

  // coupling with other physical PKs
  vapor_diffusion_ = false;
}


/* ******************************************************************
* Set defaults parameters. It should be called once
****************************************************************** */
void Richards_PK::InitializeUpwind_()
{
  darcy_flux_copy = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("darcy_flux", passwd_)));

  // Create RCP pointer to upwind flux.
  if (rel_perm_->method() == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX ||
      rel_perm_->method() == FLOW_RELATIVE_PERM_AMANZI_ARTIFICIAL_DIFFUSION ||
      rel_perm_->method() == FLOW_RELATIVE_PERM_AMANZI_MFD) {
    darcy_flux_upwind = darcy_flux_copy;
  } else if (rel_perm_->method() == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    darcy_flux_upwind = Teuchos::rcp(new CompositeVector(*darcy_flux_copy));
    rel_perm_->ComputeGravityFlux(K, gravity_, darcy_flux_upwind);
  } else {
    darcy_flux_upwind = Teuchos::rcp(new CompositeVector(*darcy_flux_copy));
    darcy_flux_upwind->PutScalar(0.0);
  }
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
  const Epetra_MultiVector& p_cell = *pressure.ViewComponent("cell");
  Epetra_MultiVector& p_face = *pressure.ViewComponent("face");

  DeriveFaceValuesFromCellValues(p_cell, p_face);

  double T1 = T_physics, T0 = T1 - dT;
  UpdateSourceBoundaryData(T0, T1, pressure);

  // saturations
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p_cell, ws);
  ws_prev = ws;
}


/* ******************************************************************
* Initial pressure is set to the pressure for fully saturated rock.
****************************************************************** */
void Richards_PK::InitializeSteadySaturated()
{ 
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "initializing with a saturated steady state..." << std::endl;
  }
  double T0 = S_->time();
  SolveFullySaturatedProblem(T0, *solution, ti_specs->solver_name_ini);
}


/* ******************************************************************
* Specific initialization of the initial pressure guess calculation.
****************************************************************** */
void Richards_PK::InitPicard(double T0)
{
  ti_igs_list_ = rp_list_->sublist("initial guess pseudo time integrator");
  ProcessSublistTimeInterval(ti_igs_list_,  ti_specs_igs_);

  ti_specs = &ti_specs_igs_;
  if (ti_specs->error_control_options == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;
  }

  InitNextTI(T0, 0.0, ti_specs_igs_);

  // calculate initial guess: cleaning is required (lipnikov@lanl.gov)
  T_physics = ti_specs_igs_.T0;
  dT = ti_specs_igs_.dT0;
  AdvanceToSteadyState_Picard(ti_specs_igs_);

  const Epetra_MultiVector& p_cell = *S_->GetFieldData("pressure")->ViewComponent("cell");
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell");

  DeriveSaturationFromPressure(p_cell, ws);
  ws_prev = ws;
}


/* ******************************************************************
* Specific initialization of a steady state time integration phase.
****************************************************************** */
void Richards_PK::InitSteadyState(double T0, double dT0)
{
  ti_sss_list_ = rp_list_->sublist("steady state time integrator");
  ProcessSublistTimeInterval(ti_sss_list_,  ti_specs_sss_);

  if (ti_specs != NULL) OutputTimeHistory(*rp_list_, ti_specs->dT_history);
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
  ti_trs_list_ = rp_list_->sublist("transient time integrator");
  ProcessSublistTimeInterval(ti_trs_list_,  ti_specs_trs_);

  if (ti_specs != NULL) OutputTimeHistory(*rp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_trs_;

  error_control_ = ti_specs->error_control_options;
  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4
    ti_specs->error_control_options = error_control_;
  }

  InitNextTI(T0, dT0, ti_specs_trs_);

  // reset some quantities
  mass_bc = 0.0;
  seepage_mass_ = 0.0;
}


/* ******************************************************************
* TBW
****************************************************************** */
void Richards_PK::InitTimeInterval()
{
  ProcessSublistTimeInterval(ti_list_, ti_specs_generic_);
 
  ti_specs_generic_.T0  = ti_list_.get<double>("start interval time", 0.0);
  ti_specs_generic_.dT0 = ti_list_.get<double>("initial time step", 1.0);

  double T0 = ti_specs_generic_.T0;
  double dT0 = ti_specs_generic_.dT0;

  dT = dT0;
  dTnext = dT0;

  if (ti_specs != NULL) OutputTimeHistory(*rp_list_, ti_specs->dT_history);
  ti_specs = &ti_specs_generic_;

  error_control_ = ti_specs->error_control_options;
  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
    ti_specs->error_control_options = error_control_;
  }

  InitNextTI(T0, dT0, ti_specs_generic_);
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Richards_PK::InitNextTI(double T0, double dT0, TI_Specs& ti_specs)
{
  ResetPKtimes(T0, dT0);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
        << vo_->color("green") << "New TI phase: " << ti_specs.ti_method_name.c_str() 
        << vo_->reset() << std::endl << std::endl;
    *vo_->os()<< "EC:" << error_control_ << " Src:" << src_sink_distribution
              << " Upwind:" << rel_perm_->method() << op_matrix_diff_->upwind()
              << " PC:\"" << ti_specs.preconditioner_name.c_str() << "\"" << std::endl;
  }

  // repeat initialization steps, mainly for old MPC.
  InitializeUpwind_();

  // set up new time integration or solver
  std::string ti_method_name(ti_specs.ti_method_name);

  if (ti_specs.ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList& bdf1_list = ti_specs.ti_list_ptr_->sublist("BDF1");
    if (! bdf1_list.isSublist("VerboseObject"))
        bdf1_list.sublist("VerboseObject") = rp_list_->sublist("VerboseObject");

    bdf1_dae = Teuchos::rcp(new BDF1_TI<CompositeVector, CompositeVectorSpace>(*this, bdf1_list, solution));
  }

  // initialize matrix and preconditioner operators
  SetAbsolutePermeabilityTensor();

  op_matrix_->Init();
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
  op_matrix_diff_->SetBCs(op_bc_);
  op_matrix_diff_->Setup(Kptr, rel_perm_->Krel(), rel_perm_->dKdP(), rho_, mu_);
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true);

  op_preconditioner_->Init();
  op_preconditioner_->SetBCs(op_bc_);
  op_preconditioner_diff_->Setup(Kptr, rel_perm_->Krel(), rel_perm_->dKdP(), rho_, mu_);
  op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->ApplyBCs(true);
  op_preconditioner_->SymbolicAssembleMatrix();

  // potentially wrap the preconditioner in a linear solver, i.e. for Newton
  if (ti_specs.solver_name != "none") {
    AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
    op_pc_solver_ = sfactory.Create(ti_specs.solver_name, *linear_operator_list_, op_preconditioner_);
  } else {
    op_pc_solver_ = op_preconditioner_;
  }
  
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    int missed_tmp = missed_bc_faces_;
    int dirichlet_tmp = dirichlet_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
    mesh_->get_comm()->SumAll(&dirichlet_tmp, &dirichlet_bc_faces_, 1);
#endif

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl << std::endl;
  }

  // Well modeling
  if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
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

  // Get a hydrostatic solution consistent with b.c.
  // Call this initialization procedure only once due to possible
  // multiple restarts of a transient time integrator.
  if (ti_specs.initialize_with_darcy) {
    SolveFullySaturatedProblem(T0, *solution, ti_specs.solver_name_ini);
    ti_specs.initialize_with_darcy = false;

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

    Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell", false);
    Epetra_MultiVector& ws_prev = *S_->GetFieldData("prev_water_saturation", passwd_)->ViewComponent("cell", false);
    DeriveSaturationFromPressure(p, ws);
    ws_prev = ws;
  }

  // initialize saturation
  Epetra_MultiVector& ws = *S_->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  DeriveSaturationFromPressure(pstate, ws);
 
  // derive mass flux (state may not have it at time 0)
  double tmp;
  darcy_flux_copy->Norm2(&tmp);
  if (tmp == 0.0) {
    op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);
  }

  // re-initialize lambda
  double T1 = T0 + dT0;
  if (flag_face && ti_specs.pressure_lambda_constraints) {
    EnforceConstraints(T1, *solution);
    // update mass flux
    op_matrix_->Init();
    op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
    op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);
  } else {
    CompositeVector& pressure = *S_->GetFieldData("pressure", passwd_);
    UpdateSourceBoundaryData(T0, T1, *solution);
    rel_perm_->Compute(pressure);

    RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
    upwind_->Compute(*darcy_flux_upwind, *solution, bc_model, bc_value,
                     *rel_perm_->Krel(), *rel_perm_->Krel(), func1);

    RelativePermeabilityUpwindFn func2 = &RelativePermeability::Derivative;
    upwind_->Compute(*darcy_flux_upwind, *solution, bc_model, bc_value,
                     *rel_perm_->dKdP(), *rel_perm_->dKdP(), func2);

    if (ti_specs.inflow_krel_correction) {
      Epetra_MultiVector& k_face = *rel_perm_->Krel()->ViewComponent("face", true);
      AmanziMesh::Entity_ID_List cells;

      for (int f = 0; f < nfaces_wghost; f++) {
        if ((bc_model[f] == Operators::OPERATOR_BC_NEUMANN || 
             bc_model[f] == Operators::OPERATOR_BC_MIXED) && bc_value[f] < 0.0) {
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

          const AmanziGeometry::Point& normal = mesh_->face_normal(f);
          double area = mesh_->face_area(f);
          double Knn = ((K[cells[0]] * normal) * normal) / (area * area);
          k_face[0][f] = std::min(1.0, -bc_value[f]  * mu_ / (Knn * rho_ * rho_ * g_));
        } 
      }    
    }
  }

  // normalize to obtain Darcy flux
  Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;

  // nonlinear solver control options
  ti_specs.num_itrs = 0;
  block_picard = 0;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    VV_PrintHeadExtrema(*solution);
  }
}


/* ******************************************************************
* This routine avoid limitations of MPC by bumping up the time step.                                          
****************************************************************** */
double Richards_PK::get_dt()
{
  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD && block_picard == 1) dT *= 1e+4;
  return dT;
}


/* ******************************************************************* 
* Performs one time step of size dT_MPC either for steady-state or 
* transient calculations.
*
* NOTE: This might require refactor for working in a more general process
*       tree. Semantic of Advance() is that it must take step of size
*       dT_MPC, otherwise return true (for failed). Steps should not be
*       taken internally in case coupling fails at a higher level.
******************************************************************* */
bool Richards_PK::Advance(double dT_MPC, double& dT_actual)
{
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


  bool fail = false;

  if (!new_mpc_driver){
    if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
      while (bdf1_dae->TimeStep(dT, dTnext, solution)) {
	dT = dTnext;
      }
      // --etc this is a bug in general -- should not commit the solution unless
      //   the step passes all other PKs
      bdf1_dae->CommitSolution(dT, solution);
      T_physics = bdf1_dae->time();
    }

    dT_actual = dT;
  }
  else {
    if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1){
      fail = bdf1_dae->TimeStep(dT, dTnext, solution);
      if (fail){
	dT = dTnext;
	return fail;
      }
      bdf1_dae->CommitSolution(dT, solution);
      T_physics = bdf1_dae->time();
    }
  }


  // tell the caller what time step we actually took
  //dT_actual = dT;
  
  dt_tuple times(time, dT);
  ti_specs->dT_history.push_back(times);

  ti_specs->num_itrs++;


 if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
   VV_ReportWaterBalance(S_.ptr());
  }
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    VV_ReportSeepageOutflow(S_.ptr());
  }


  dT = dTnext;

  //exit(0);
  
  return fail;
}


/* ******************************************************************
* Transfer part of the internal data needed by transport to the 
* flow state FS_MPC. MPC may request to populate the original FS.
* The consistency condition is improved by adjusting saturation while
* preserving its LED property.
****************************************************************** */
void Richards_PK::CommitState(double dt, const Teuchos::Ptr<State>& S)
{
  //UpdateAuxilliaryData();

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
  op_matrix_diff_->UpdateFlux(*solution, darcy_flux);

  Epetra_MultiVector& flux = *darcy_flux.ViewComponent("face", true);
  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;
  *darcy_flux_copy->ViewComponent("face", true) = flux;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  // update mass balance
  // ImproveAlgebraicConsistency(ws_prev, ws);
 
  dT = dTnext;
}


/* ******************************************************************
 * * A wrapper for updating boundary conditions.
 * ****************************************************************** */
void Richards_PK::UpdateSourceBoundaryData(double T0, double T1, const CompositeVector& u)
{
  if (src_sink != NULL) {
    if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(T0, T1, Kxy->Values());
    } else {
      src_sink->ComputeDistribute(T0, T1, NULL);
    }
  }

  bc_pressure->Compute(T1);
  bc_flux->Compute(T1);
  bc_seepage->Compute(T1);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(T1);
  else
    bc_head->ComputeShift(T1, shift_water_table_->Values());

  ComputeBCs(u);
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
  S_->GetFieldData("water_saturation", passwd_)->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);
  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell", false);

  WhetStone::MFD3D_Diffusion mfd3d(mesh_);
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
      int c2 = mfd3d.cell_get_face_adj_cell(c, f);
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



/* ******************************************************************
* 
****************************************************************** */
double Richards_PK::BoundaryFaceValue(int f, const CompositeVector& u)
{
  const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
  double face_value;

  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    face_value = u_face[0][f];
  } else {
    Epetra_MultiVector& k_face = *rel_perm_->Krel()->ViewComponent("face");
    Epetra_MultiVector& dk_face = *rel_perm_->dKdP()->ViewComponent("face");

    std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_->WRM();
    const Epetra_IntVector& map_c2mb = rel_perm_->map_c2mb();

    int c = BoundaryFaceGetCell(f);
    face_value = DeriveBoundaryFaceValue(f, u, WRM[map_c2mb[c]]);
    // k_face[0][f] = WRM[map_c2mb[c]]->k_relative (atm_pressure_ - face_val);
    // dk_face[0][f] = -WRM[map_c2mb[c]]->dKdPc (atm_pressure_ - face_val);
  }
  return face_value;
}

void  Richards_PK::CalculateDiagnostics(const Teuchos::Ptr<State>& S){
  UpdateAuxilliaryData();
}

}  // namespace Flow
}  // namespace Amanzi

