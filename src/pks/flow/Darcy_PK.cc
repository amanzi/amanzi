/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <vector>

// TPLs
#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Amanzi
#include "CommonDefs.hh"
#include "EvaluatorPrimary.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "MeshAlgorithms.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_Utils.hh"
#include "TimestepControllerFactory.hh"
#include "Tensor.hh"

// Amanzi::Flow
#include "ApertureDarcyEvaluator.hh"
#include "Darcy_PK.hh"
#include "FlowDefs.hh"
#include "ModelEvaluator.hh"
#include "SpecificStorage.hh"
#include "StateHelpers.hh"
#include "VolumetricFlowRateEvaluator.hh"
#include "WaterStorageDarcy.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
Darcy_PK::Darcy_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln), Flow_PK(pk_tree, glist, S, soln), soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator", true);

  // domain and primary evaluators
  domain_ = fp_list_->template get<std::string>("domain name", "domain");
  pressure_key_ = Keys::getKey(domain_, "pressure");
  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");

  AddDefaultPrimaryEvaluator(S_, pressure_key_);
  AddDefaultPrimaryEvaluator(S_, mol_flowrate_key_);
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
Darcy_PK::Darcy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const std::string& pk_list_name,
                   Teuchos::RCP<State> S,
                   const Teuchos::RCP<TreeVector>& soln)
  : Flow_PK(), soln_(soln)
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_list_name, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator", true);

  // domain and primary evaluators
  domain_ = fp_list_->template get<std::string>("domain name", "domain");
  pressure_key_ = Keys::getKey(domain_, "pressure");
  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");

  AddDefaultPrimaryEvaluator(S_, pressure_key_);
  AddDefaultPrimaryEvaluator(S_, mol_flowrate_key_);
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void
Darcy_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);
  dim = mesh_->getSpaceDimension();

  // generate keys used by Darcy PK only
  specific_yield_key_ = Keys::getKey(domain_, "specific_yield");
  specific_storage_key_ = Keys::getKey(domain_, "specific_storage");

  ref_aperture_key_ = Keys::getKey(domain_, "ref_aperture");
  ref_pressure_key_ = Keys::getKey(domain_, "ref_pressure");
  compliance_key_ = Keys::getKey(domain_, "compliance");

  // optional keys
  pressure_head_key_ = Keys::getKey(domain_, "pressure_head");

  // set up the base class
  Flow_PK::Setup();

  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(fp_list_, "physical models and assumptions");
  external_aperture_ = physical_models->get<bool>("external aperture", false);
  std::string mu_model = physical_models->get<std::string>("viscosity model", "constant viscosity");
  use_bulk_modulus_ = physical_models->get<bool>("use bulk modulus", false);
  if (mu_model != "constant viscosity") {
    Errors::Message msg;
    msg << "Darcy PK supports only constant viscosity model.";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::sublist(fp_list_, "operators", true);
  Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::sublist(list1, "diffusion operator", true);
  Teuchos::RCP<Teuchos::ParameterList> list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  // primary field: pressure
  if (!S_->HasRecord(pressure_key_)) {
    std::vector<std::string> names;
    std::vector<AmanziMesh::Entity_kind> locations;
    std::vector<int> ndofs;

    names.push_back("cell");
    locations.push_back(AmanziMesh::Entity_kind::CELL);
    ndofs.push_back(1);

    if (name != "fv: default" && name != "nlfv: default") {
      names.push_back("face");
      locations.push_back(AmanziMesh::Entity_kind::FACE);
      ndofs.push_back(1);
    }

    S_->Require<CV_t, CVS_t>(pressure_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  // require additional fields and evaluators
  // Many fields/evaluators have a simple struncture. They are ghosted
  //   cell-based fields. We use a helper function that reruires a field
  //   and returns a parameter list populated with standard values.
  // -- water storage
  if (!S_->HasRecord(water_storage_key_)) {
    auto elist = RequireFieldForEvaluator(*S_, water_storage_key_);
    S_->GetRecordW(water_storage_key_, water_storage_key_).set_io_vis(false);

    elist.set<std::string>("pressure key", pressure_key_)
      .set<std::string>("specific storage key", specific_storage_key_);
    if (flow_on_manifold_) elist.set<std::string>("aperture key", aperture_key_);

    S_->RequireDerivative<CV_t, CVS_t>(
        water_storage_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, water_storage_key_)
      .SetGhosted();

    auto eval = Teuchos::rcp(new WaterStorageDarcy(elist));
    S_->SetEvaluator(water_storage_key_, Tags::DEFAULT, eval);
  }

  // -- water storage from the previous time step
  if (!S_->HasRecord(prev_water_storage_key_)) {
    S_->Require<CV_t, CVS_t>(prev_water_storage_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->GetRecordW(prev_water_storage_key_, passwd_).set_io_vis(false);
  }

  // -- specific storage
  if (!S_->HasRecord(specific_storage_key_)) {
    RequireFieldForEvaluator(*S_, specific_storage_key_);
    S_->RequireEvaluator(specific_storage_key_, Tags::DEFAULT);
  }

  if (!S_->HasRecord(specific_yield_key_)) {
    S_->Require<CV_t, CVS_t>(specific_yield_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- saturation
  if (!S_->HasRecord(saturation_liquid_key_)) {
    RequireFieldForEvaluator(*S_, saturation_liquid_key_);
    S_->RequireEvaluator(saturation_liquid_key_, Tags::DEFAULT);
  }

  if (!S_->HasRecord(prev_saturation_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(prev_saturation_liquid_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- porosity
  if (!S_->HasRecord(porosity_key_)) {
    RequireFieldForEvaluator(*S_, porosity_key_);
    S_->RequireEvaluator(porosity_key_, Tags::DEFAULT);
  }

  // -- viscosity
  if (!S_->HasRecord("const_fluid_viscosity")) {
    S_->Require<double>("const_fluid_viscosity", Tags::DEFAULT, "state");
  }

  // -- molar and volumetric flow rates
  double rho = S_->ICList().sublist("const_fluid_density").get<double>("value");
  double molar_rho = rho / CommonDefs::MOLAR_MASS_H2O;
  Setup_FlowRates_(true, molar_rho);

  // -- fracture dynamics
  if (flow_on_manifold_) {
    S_->Require<CV_t, CVS_t>(compliance_key_, Tags::DEFAULT, compliance_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(compliance_key_, Tags::DEFAULT);

    if (external_aperture_) {
      RequireFieldForEvaluator(*S_, ref_aperture_key_);
      S_->GetRecordW(ref_aperture_key_, ref_aperture_key_).set_io_vis(false);
      S_->RequireEvaluator(ref_aperture_key_, Tags::DEFAULT);

      S_->Require<CV_t, CVS_t>(ref_pressure_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      S_->GetRecordW(ref_pressure_key_, passwd_).set_io_vis(false);
      S_->GetRecordW(ref_pressure_key_, passwd_).set_io_checkpoint(false);

      S_->RequireEvaluator(aperture_key_, Tags::DEFAULT);
    } else {
      S_->RequireEvaluator(aperture_key_, Tags::DEFAULT);
    }
  } else if (use_bulk_modulus_) {
    S_->Require<CV_t, CVS_t>(bulk_modulus_key_, Tags::DEFAULT)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(bulk_modulus_key_, Tags::DEFAULT);
  }

  // Require additional components for the existing fields
  Teuchos::ParameterList abs_perm = fp_list_->sublist("absolute permeability");
  coordinate_system_ = abs_perm.get<std::string>("coordinate system", "cartesian");
  int noff = abs_perm.get<int>("off-diagonal components", 0);

  if (noff > 0) {
    CompositeVectorSpace& cvs =
      S_->Require<CV_t, CVS_t>(permeability_key_, Tags::DEFAULT, permeability_key_);
    cvs.SetOwned(false);
    cvs.AddComponent("offd", AmanziMesh::Entity_kind::CELL, noff)->SetOwned(true);
  }

  // require optional fields
  if (fp_list_->isParameter("optional fields")) {
    std::vector<std::string> fields =
      fp_list_->get<Teuchos::Array<std::string>>("optional fields").toVector();
    for (auto it = fields.begin(); it != fields.end(); ++it) {
      Key optional_key = Keys::getKey(domain_, *it);
      if (!S_->HasRecord(optional_key)) {
        S_->Require<CV_t, CVS_t>(optional_key, Tags::DEFAULT, passwd_)
          .SetMesh(mesh_)
          ->SetGhosted(true)
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    }
  }

  // to force other PKs to use density, we define it here
  if (S_->HasEvaluator(mol_density_liquid_key_, Tags::DEFAULT))
    AddDefaultIndependentEvaluator(S_, mol_density_liquid_key_, Tags::DEFAULT, molar_rho);

  if (!S_->HasRecord(mass_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT, mass_density_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    AddDefaultIndependentEvaluator(S_, mass_density_liquid_key_, Tags::DEFAULT, rho);
  }

  // -- hyadraulic head and full Darcy velocity
  Setup_LocalFields_();

  // save frequently used evaluators
  pressure_eval_ = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(pressure_key_, Tags::DEFAULT));

  // set units
  S_->GetRecordSetW(pressure_key_).set_units("Pa");
  S_->GetRecordSetW(specific_storage_key_).set_units("m");
  S_->GetRecordSetW(hydraulic_head_key_).set_units("m");
  S_->GetRecordSetW(porosity_key_).set_units("-");
  S_->GetRecordSetW(saturation_liquid_key_).set_units("-");
  S_->GetRecordSetW(hydraulic_head_key_).set_units("m");
}


/* ******************************************************************
* Extract information from parameter list and initialize data.
****************************************************************** */
void
Darcy_PK::Initialize()
{
  // Initialize miscalleneous defaults.
  // -- others
  initialize_with_darcy_ = true;
  num_itrs_ = 0;

  // Create verbosity object to print out initialization statisticsr.,
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = fp_list_->sublist("verbose object");

  std::string ioname = "DarcyPK";
  if (domain_ != "domain") ioname += "-" + domain_;
  vo_ = Teuchos::rcp(new VerboseObject(ioname, vlist));

  // Initilize various base class data.
  Flow_PK::Initialize();

  // -- times
  double t_ini = S_->get_time();
  dt_next_ = dt_;
  dt_desirable_ = dt_; // The minimum desirable time step from now on.
  dt_history_.clear();

  // Initialize local fields and evaluators.
  InitializeFields_();
  UpdateLocalFields_(S_.ptr());

  // Create solution and auxiliary data for time history.
  solution = S_->GetPtrW<CV_t>(pressure_key_, Tags::DEFAULT, passwd_);
  soln_->SetData(solution);

  const Epetra_BlockMap& cmap = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));

  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  if (ti_method_name == "BDF1") {
    Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE; // usually 1e-4;

    // time step controller
    TimestepControllerFactory<Epetra_MultiVector> fac;
    ts_control_ = fac.Create(bdf1_list, pdot_cells, pdot_cells_prev);
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: BDF1 time integration list is missing..." << std::endl;
  }

  // Initialize specific yield
  specific_yield_copy_ = Teuchos::null;
  UpdateSpecificYield_();

  // Initialize lambdas. It may be used by boundary conditions.
  auto& pressure = S_->GetW<CompositeVector>(pressure_key_, Tags::DEFAULT, passwd_);

  if (ti_list_->isSublist("pressure-lambda constraints") && solution->HasComponent("face")) {
    std::string method =
      ti_list_->sublist("pressure-lambda constraints").get<std::string>("method");
    if (method == "projection") {
      Epetra_MultiVector& p = *solution->ViewComponent("cell");
      Epetra_MultiVector& lambda = *solution->ViewComponent("face");

      DeriveFaceValuesFromCellValues(p, lambda);
      S_->GetRecordW(pressure_key_, Tags::DEFAULT, passwd_).set_initialized(true);
    }
  }

  // Create and initialize boundary conditions and source terms.
  flux_units_ = 1.0;
  InitializeBCsSources_(*fp_list_);
  UpdateSourceBoundaryData(t_ini, t_ini, pressure);

  // Initialize diffusion operator and solver.
  // -- instead of scaling K, we scale the elemental mass matrices
  double mu = S_->Get<double>("const_fluid_viscosity");
  Teuchos::ParameterList& oplist =
    fp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");
  if (flow_on_manifold_) oplist.set<std::string>("nonlinear coefficient", "standard: cell");
  if (coupled_to_matrix_ || flow_on_manifold_) {
    if (!oplist.isParameter("use manifold flux")) oplist.set<bool>("use manifold flux", true);
  }

  Operators::PDE_DiffusionFactory opfactory(oplist, mesh_);
  opfactory.SetConstantGravitationalTerm(gravity_, rho_);

  if (!flow_on_manifold_) {
    SetAbsolutePermeabilityTensor();
    Teuchos::RCP<std::vector<WhetStone::Tensor>> Kptr = Teuchos::rcpFromRef(K);
    opfactory.SetVariableTensorCoefficient(Kptr);
    opfactory.SetConstantScalarCoefficient(rho_ / mu);
  } else {
    WhetStone::Tensor Ktmp(dim, 1);
    Ktmp(0, 0) = rho_ / mu;
    opfactory.SetConstantTensorCoefficient(Ktmp);

    S_->GetEvaluator(permeability_eff_key_).Update(*S_, permeability_eff_key_);
    auto kptr = S_->GetPtr<CV_t>(permeability_eff_key_, Tags::DEFAULT);
    opfactory.SetVariableScalarCoefficient(kptr);
  }

  op_diff_ = opfactory.Create();
  op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_diff_->SetBCs(op_bc_, op_bc_);
  op_ = op_diff_->global_operator();

  // -- accumulation operator
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op_));
  op_->CreateCheckPoint();

  // -- generic linear solver.
  AMANZI_ASSERT(ti_list_->isParameter("linear solver"));
  solver_name_ = ti_list_->get<std::string>("linear solver");

  // -- preconditioner. There is no need to enhance it for Darcy
  AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  std::string name = ti_list_->get<std::string>("preconditioner");
  op_->set_inverse_parameters(
    name, *preconditioner_list_, solver_name_, *linear_operator_list_, true);
  op_->InitializeInverse();

  // Optional step: calculate hydrostatic solution consistent with BCs.
  // We have to do it only once per time period.
  bool init_darcy(false);
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_) {
    initialize_with_darcy_ = false;
    bool wells_on = ti_list_->sublist("initialization").get<bool>("active wells", false);
    SolveFullySaturatedProblem(*solution, wells_on);
    init_darcy = true;
  }

  // set up operators for evaluators
  auto eval = S_->GetEvaluatorPtr(vol_flowrate_key_, Tags::DEFAULT);
  Teuchos::rcp_dynamic_cast<VolumetricFlowRateEvaluator>(eval)->set_bc(op_bc_);

  // Verbose output of initialization statistics.
  InitializeStatistics_(init_darcy);
}


/* ****************************************************************
* This completes initialization of common fields that were not
* initialized by the state.
**************************************************************** */
void
Darcy_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  InitializeCVField(S_, *vo_, saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_key_, 1.0);
  InitializeCVField(S_, *vo_, prev_saturation_liquid_key_, Tags::DEFAULT, passwd_, 1.0);

  if (flow_on_manifold_)
    InitializeCVField(S_, *vo_, compliance_key_, Tags::DEFAULT, compliance_key_, 0.0);

  if (flow_on_manifold_ && external_aperture_)
    InitializeCVFieldFromCVField(S_, *vo_, ref_pressure_key_, pressure_key_, passwd_);

  InitializeCVFieldFromCVField(S_, *vo_, prev_water_storage_key_, water_storage_key_, passwd_);
}


/* ******************************************************************
* Print the header for new time period.
****************************************************************** */
void
Darcy_PK::InitializeStatistics_(bool init_darcy)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    double mu = S_->Get<double>("const_fluid_viscosity");
    std::string ti_name = ti_list_->get<std::string>("time integration method", "none");
    std::string ts_name =
      (ti_name == "BDF1") ?
        ti_list_->sublist(ti_name).get<std::string>("timestep controller type") :
        "timestep controller is not defined";
    std::string pc_name = ti_list_->get<std::string>("preconditioner");

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nTI:\"" << ti_name.c_str() << "\""
               << " dt:" << ts_name << " LS:\"" << solver_name_.c_str() << "\""
               << " PC:\"" << pc_name.c_str() << "\"" << std::endl
               << "matrix: " << op_->PrintDiagnostics() << std::endl;
    *vo_->os() << "constant viscosity model, mu=" << mu << std::endl;

    if (init_darcy) {
      *vo_->os() << "initial pressure guess: \"from saturated solver\"\n" << std::endl;
    } else {
      *vo_->os() << "initial pressure guess: \"from State\"\n" << std::endl;
    }

    int missed_tmp = missed_bc_faces_;
    int dirichlet_tmp = dirichlet_bc_faces_;
#ifdef HAVE_MPI
    mesh_->getComm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
    mesh_->getComm()->SumAll(&dirichlet_tmp, &dirichlet_bc_faces_, 1);
#endif

    *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl
               << std::endl;

    VV_PrintHeadExtrema(*solution);
    VV_PrintSourceExtrema();

    *vo_->os() << vo_->color("green") << "Initialization of PK is complete, T=" << S_->get_time()
               << " dT=" << get_dt() << vo_->reset() << std::endl
               << std::endl;
  }

  if (dirichlet_bc_faces_ == 0 && domain_ == "domain" && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }
}


/* *******************************************************************
* Performs one time step from t_old to t_new. The boundary conditions
* are calculated only once, during the initialization step.
******************************************************************* */
bool
Darcy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;
  double dt_MPC(dt_);

  // refresh data
  UpdateSourceBoundaryData(t_old, t_new, *solution);

  // add accumulation term: specific storage
  op_->RestoreCheckPoint();

  const auto& ws_prev = S_->Get<CV_t>(prev_water_storage_key_, Tags::DEFAULT);
  S_->GetEvaluator(water_storage_key_).UpdateDerivative(*S_, passwd_, pressure_key_, Tags::DEFAULT);
  const auto& dws_dp =
    S_->GetDerivative<CV_t>(water_storage_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT);
  op_acc_->AddAccumulationRhs(dws_dp, ws_prev, 1.0 / dt_, "cell", true);
  // op_acc_->AddAccumulationDelta(*solution, dws_dp, dws_dp, dt_, "cell");

  // add accumulation term: specific yield
  double factor = 1.0 / (g_ * dt_);
  CompositeVector sy_g(*specific_yield_copy_);
  sy_g.Scale(factor);

  if (flow_on_manifold_) {
    S_->GetEvaluator(aperture_key_).Update(*S_, aperture_key_);
    const auto& aperture = S_->Get<CV_t>(aperture_key_, Tags::DEFAULT);
    sy_g.Multiply(1.0, sy_g, aperture, 0.0);
  }

  op_acc_->AddAccumulationDeltaNoVolume(*solution, sy_g, "cell");

  // Peaceman model
  if (S_->HasRecord("well_index")) {
    const auto& wi = S_->Get<CV_t>("well_index");
    op_acc_->AddAccumulationTerm(wi, "cell");
  }

  // add diffusion matrices
  if (flow_on_manifold_) {
    S_->GetEvaluator(permeability_eff_key_).Update(*S_, "flow");
    op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  }
  op_diff_->ApplyBCs(true, true, true);

  CompositeVector& rhs = *op_->rhs();
  AddSourceTerms(rhs, dt_);

  op_->ComputeInverse();

  // save pressure at time t^n.
  std::string dt_control = ti_list_->sublist("BDF1").get<std::string>("timestep controller type");
  Teuchos::RCP<Epetra_MultiVector> p_old;
  if (dt_control == "adaptive") {
    p_old = Teuchos::rcp(new Epetra_MultiVector(*solution->ViewComponent("cell")));
  }

  op_->ApplyInverse(rhs, *solution);
  pressure_eval_->SetChanged();

  // statistics
  num_itrs_++;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    double pnorm;
    solution->Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver_name_ << "): ||p,lambda||=" << pnorm
               << "  itrs=" << op_->num_itrs() << std::endl;
    VV_PrintHeadExtrema(*solution);
  }

  // calculate time derivative and 2nd-order solution approximation
  if (dt_control == "adaptive") {
    Epetra_MultiVector& p_new = *solution->ViewComponent("cell"); // pressure at t^{n+1}

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = (p_new[0][c] - (*p_old)[0][c]) / dt_;
      p_new[0][c] = (*p_old)[0][c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dt_ / 2;
    }
  }

  // estimate time multiplier
  dt_desirable_ = ts_control_->get_timestep(dt_MPC, 1);

  // Darcy_PK always takes the suggested time step and cannot fail
  dt_tuple times(t_new, dt_MPC);
  dt_history_.push_back(times);

  return false;
}


/* ******************************************************************
* Transfer data from the external flow state FS_MPC. MPC may request
* to populate the original state FS.
****************************************************************** */
void
Darcy_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  ComputeMolarFlowRate_(true);
  S_->GetEvaluator(vol_flowrate_key_).Update(*S_, passwd_);

  S_->GetEvaluator(water_storage_key_).Update(*S_, "flow");
  S_->GetW<CV_t>(prev_water_storage_key_, Tags::DEFAULT, passwd_) =
    S_->Get<CV_t>(water_storage_key_, Tags::DEFAULT);

  if (coupled_to_matrix_ || flow_on_manifold_) VV_FractureConservationLaw();

  // update time derivative
  *pdot_cells_prev = *pdot_cells;
}


/* ******************************************************************
* Add area/length factor to specific yield.
****************************************************************** */
void
Darcy_PK::UpdateSpecificYield_()
{
  specific_yield_copy_ =
    Teuchos::rcp(new CompositeVector(S_->Get<CV_t>(specific_yield_key_), true));

  // do we have non-zero specific yield?
  double tmp;
  specific_yield_copy_->Norm2(&tmp);
  if (tmp == 0.0) return;

  // populate ghost cells
  specific_yield_copy_->ScatterMasterToGhosted();
  Epetra_MultiVector& specific_yield = *specific_yield_copy_->ViewComponent("cell", true);

  int negative_yield = 0;

  for (int c = 0; c < ncells_owned; c++) {
    if (specific_yield[0][c] > 0.0) {
      auto [faces, dirs] = mesh_->getCellFacesAndDirections(c);
      auto adjcells =
        AmanziMesh::getCellFaceAdjacentCells(*mesh_, c, AmanziMesh::Parallel_kind::OWNED);

      double area = 0.0;
      int nfaces = faces.size();
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        int c2 = adjcells[f];
        if (c2 >= 0) {
          if (specific_yield[0][c2] <= 0.0) // cell in the fully saturated layer
            area -= (mesh_->getFaceNormal(f))[dim - 1] * dirs[n];
        }
      }
      specific_yield[0][c] *= area;
      if (area <= 0.0) negative_yield++;
    }
  }

#ifdef HAVE_MPI
  int negative_yield_tmp = negative_yield;
  mesh_->getComm()->MaxAll(&negative_yield_tmp, &negative_yield, 1);
#endif
  if (negative_yield > 0) {
    Errors::Message msg;
    msg << "Flow PK: configuration of the yield region leads to negative yield interfaces.";
    Exceptions::amanzi_throw(msg);
  }
}


} // namespace Flow
} // namespace Amanzi
