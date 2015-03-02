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
#include "independent_variable_field_evaluator_fromfunction.hh"
#include "Mesh.hh"
#include "mfd3d_diffusion.hh"
#include "OperatorDiffusionFactory.hh"
#include "Point.hh"
#include "primary_variable_field_evaluator.hh"
#include "UpwindFactory.hh"
#include "XMLParameterListWriter.hh"

#include "darcy_velocity_evaluator.hh"
#include "Flow_BC_Factory.hh"
#include "Richards_PK.hh"
#include "WRMEvaluator.hh"

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
  
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);
  ti_list_ = Teuchos::sublist(rp_list_, "time integrator");
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
    pressure_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("pressure", pressure_eval_);
  }
  if (!S_->HasField("hydraulic_head")) {
    S_->RequireField("hydraulic_head", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);
  }

  if (!S_->HasField("prev_saturation_liquid")) {
    S_->RequireField("prev_saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetField("prev_saturation_liquid", passwd_)->set_io_vis(false);
  }

  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "darcy_flux");
    darcy_flux_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("darcy_flux", darcy_flux_eval_);
  }

  if (!S_->HasField("darcy_velocity")) {
    S_->RequireField("darcy_velocity", "darcy_velocity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist;
    Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S_->SetFieldEvaluator("darcy_velocity", eval);
  }

  // Additional structutors for evaluators to work properly.
  Teuchos::RCP<Teuchos::ParameterList>
      wrm_list = Teuchos::sublist(rp_list_, "water retention models", true);
  wrm_ = CreateWRMPartition(mesh_, wrm_list);

  // requesting evaluators
  if (!S_->HasField("saturation_liquid")) {
    S_->RequireField("saturation_liquid", "saturation_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.sublist("VerboseObject").set<std::string>("Verbosity Level", "extreme");
    Teuchos::RCP<WRMEvaluator> eval = Teuchos::rcp(new WRMEvaluator(elist, wrm_));
    S_->SetFieldEvaluator("saturation_liquid", eval);
  }

  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", "porosity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireFieldEvaluator("porosity");
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
  // solution = Teuchos::rcp(new CompositeVector(op_matrix_->DomainMap()));
  solution = S_->GetFieldData("pressure", passwd_);
  
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
  initialize_with_darcy_ = true;
  num_itrs_ = 0;

  // coupling with other physical PKs
  vapor_diffusion_ = false;
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
  pressure_eval_->SetFieldAsChanged(S_.ptr());
  S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
  Epetra_MultiVector& s_l = *S_->GetFieldData("saturation_liquid", "saturation_liquid")->ViewComponent("cell");

  Epetra_MultiVector& s_l_prev = *S_->GetFieldData("prev_saturation_liquid", passwd_)->ViewComponent("cell");
  s_l_prev = s_l;
}


/* ******************************************************************
* This is long but simple subroutine. It goes through time integrator
* list and initializes various objects created during setup step.
****************************************************************** */
void Richards_PK::InitTimeInterval()
{
  // times
  double T0 = ti_list_->get<double>("start interval time", 0.0);
  double dT0 = ti_list_->get<double>("initial time step", 1.0);
  ResetPKtimes(T0, dT0);

  dT = dT0;
  dTnext = dT0;

  // error control options
  ASSERT(ti_list_->isParameter("error control options"));

  error_control_ = 0;
  std::vector<std::string> options;
  options = ti_list_->get<Teuchos::Array<std::string> >("error control options").toVector();

  for (int i=0; i < options.size(); i++) {
    if (options[i] == "pressure") {
      error_control_ += FLOW_TI_ERROR_CONTROL_PRESSURE;
    } else if (options[i] == "saturation") {
      error_control_ += FLOW_TI_ERROR_CONTROL_SATURATION;
    } else if (options[i] == "residual") {
      error_control_ += FLOW_TI_ERROR_CONTROL_RESIDUAL;
    }
  }

  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
  }

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  if (! bdf1_list.isSublist("VerboseObject"))
      bdf1_list.sublist("VerboseObject") = rp_list_->sublist("VerboseObject");

  bdf1_dae = Teuchos::rcp(new BDF1_TI<CompositeVector, CompositeVectorSpace>(*this, bdf1_list, solution));

  // complete other steps
  // repeat upwind initialization, mainly for old MPC
  InitializeUpwind_();

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

  // preconditioner and optional linear solver
  ASSERT(ti_list_->isParameter("preconditioner"));
  preconditioner_name_ = ti_list_->get<std::string>("preconditioner");
  ASSERT(preconditioner_list_->isSublist(preconditioner_name_));
  
  ASSERT(ti_list_->isParameter("linear solver"));
  solver_name_ = ti_list_->get<std::string>("linear solver");

  if (solver_name_ != "none") {
    AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
    op_pc_solver_ = sfactory.Create(solver_name_, *linear_operator_list_, op_preconditioner_);
  } else {
    op_pc_solver_ = op_preconditioner_;
  }
  
  // initialize well modeling
  if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    CalculatePermeabilityFactorInWell();
  }

  // Optional step: calculate hydrostatic solution consistent with BCs
  // and clip it as requested. We have to do it only once per time period.
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_) {
    initialize_with_darcy_ = false;
    Teuchos::ParameterList& ini_list = ti_list_->sublist("initialization");
 
    std::string solver_name_ini = ini_list.get<std::string>("method", "none");
    if (solver_name_ini == "saturated solver") {
      std::string name = ini_list.get<std::string>("linear solver");
      SolveFullySaturatedProblem(T0, *solution, name);

      bool clip(false);
      double clip_saturation = ini_list.get<double>("clipping saturation value", -1.0);
      if (clip_saturation > 0.0) {
        double pmin = FLOW_PRESSURE_ATMOSPHERIC;
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(pmin, clip_saturation, p);
        clip = true;
      }

      double clip_pressure = ini_list.get<double>("clipping pressure value", -1e+10);
      if (clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(clip_pressure, p);
        clip = true;
      }

      if (clip && solution->HasComponent("face")) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        Epetra_MultiVector& lambda = *solution->ViewComponent("face", true);
        DeriveFaceValuesFromCellValues(p, lambda);
      }
    }
    else if (solver_name_ini == "picard") {
      AdvanceToSteadyState_Picard(ti_list_->sublist("initialization").sublist("picard parameters"));
    }
    pressure_eval_->SetFieldAsChanged(S_.ptr());

    // initialization is usually doen at time 0, so we need to update other
    // fields such as prev_saturation_liquid
    S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
    Epetra_MultiVector& s_l = *S_->GetFieldData("saturation_liquid", "saturation_liquid")->ViewComponent("cell");

    Epetra_MultiVector& s_l_prev = *S_->GetFieldData("prev_saturation_liquid", passwd_)->ViewComponent("cell");
    s_l_prev = s_l;
  }

  // Trigger update of secondary fields depending on the primary pressure.
  pressure_eval_->SetFieldAsChanged(S_.ptr());

  // derive mass flux (state may not have it at time 0)
  double tmp;
  darcy_flux_copy->Norm2(&tmp);
  if (tmp == 0.0) {
    op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);

    Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
    for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;
  }

  // subspace entering: re-inititime lambdas.
  double T1 = T0 + dT0;
  if (ti_list_->isSublist("pressure-lambda constraints") && solution->HasComponent("face")) {
    solver_name_constraint_ = ti_list_->sublist("pressure-lambda constraints").get<std::string>("linear solver");

    EnforceConstraints(T1, *solution);
    pressure_eval_->SetFieldAsChanged(S_.ptr());

    // update mass flux
    op_matrix_->Init();
    op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
    op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);

    // normalize to Darcy flux, m/s
    Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face", true);
    for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;

    InitializeUpwind_();
  }

  // verbose output
  // print the header for new time period
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
        << vo_->color("green") << "Initalization of TI period is complete: " << ti_method_name.c_str() 
        << vo_->reset() << std::endl << std::endl;
    *vo_->os()<< "EC:" << error_control_ << " Src:" << src_sink_distribution
              << " Upwind:" << rel_perm_->method() << op_matrix_diff_->upwind()
              << " PC:\"" << preconditioner_name_.c_str() << "\"" << std::endl;

    int missed_tmp = missed_bc_faces_;
    int dirichlet_tmp = dirichlet_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
    mesh_->get_comm()->SumAll(&dirichlet_tmp, &dirichlet_bc_faces_, 1);
#endif

    *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl << std::endl;

    VV_PrintHeadExtrema(*solution);
  }
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

  // initialization
  time = T_physics;
  if (num_itrs_ == 0) {
    Teuchos::RCP<CompositeVector> udot = Teuchos::rcp(new CompositeVector(*solution));
    udot->PutScalar(0.0);
    bdf1_dae->SetInitialState(time, solution, udot);

    UpdatePreconditioner(time, solution, dT);
    num_itrs_++;
  }

  bool fail = false;
  fail = bdf1_dae->TimeStep(dT, dTnext, solution);
  if (fail) {
    dT = dTnext;
    return fail;
  }
  bdf1_dae->CommitSolution(dT, solution);
  T_physics = bdf1_dae->time();
  pressure_eval_->SetFieldAsChanged(S_.ptr());

  // tell the caller what time step we actually took
  // dT_actual = dT;
  
  dt_tuple times(time, dT);
  dT_history_.push_back(times);
  num_itrs_++;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    VV_ReportWaterBalance(S_.ptr());
  }
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    VV_ReportSeepageOutflow(S_.ptr());
  }

  dT = dTnext;
  
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
  // ws -> ws_prev
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell");
  Epetra_MultiVector& s_l_prev = *S->GetFieldData("prev_saturation_liquid", passwd_)->ViewComponent("cell");
  s_l_prev = s_l;

  // calculate new water saturation
  S->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S.ptr(), "flow");

  // calculate Darcy flux as diffusive part + advective part.
  CompositeVector& darcy_flux = *S->GetFieldData("darcy_flux", passwd_);
  op_matrix_diff_->UpdateFlux(*solution, darcy_flux);

  Epetra_MultiVector& flux = *darcy_flux.ViewComponent("face", true);
  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;
  *darcy_flux_copy->ViewComponent("face", true) = flux;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

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

    std::vector<Teuchos::RCP<WRM> >& WRM = rel_perm_->wrm();
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

