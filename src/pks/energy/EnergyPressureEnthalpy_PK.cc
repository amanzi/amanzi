/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK

  Process kernel for energy equation in the pressure-enthalpy variables.
*/

// Amanzi
#include "CommonDefs.hh"
#include "IAPWS97.hh"
#include "LScheme_Helpers.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "UpwindFactory.hh"

// Amanzi::Energy
#include "ApertureModelEvaluator.hh"
#include "EnergyPressureEnthalpy_PK.hh"
#include "EnthalpyEvaluator.hh"
#include "EvaluatorMultiplicativeReciprocal.hh"
#include "IEMEvaluator.hh"
#include "StateArchive.hh"
#include "ThermodynamicStateEvaluators.hh"
#include "TotalEnergyEvaluatorPH.hh"

namespace Amanzi {
namespace Energy {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

static constexpr double cfactor = 1000.0 * CommonDefs::MOLAR_MASS_H2O;

/* ******************************************************************
* Default constructor.
****************************************************************** */
EnergyPressureEnthalpy_PK::EnergyPressureEnthalpy_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& glist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln),
    Energy_PK(),
    glist_(glist), 
    soln_(soln),
    passwd_("")
{
  auto pk_list = Teuchos::sublist(glist, "PKs", true);
  ep_list_ = Teuchos::sublist(pk_list, name_, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  ti_list_ = Teuchos::sublist(ep_list_, "time integrator");

  // domain name
  domain_ = ep_list_->get<std::string>("domain name", "domain");
  enthalpy_key_ = Keys::getKey(domain_, "enthalpy");
  AddDefaultPrimaryEvaluator(S_, enthalpy_key_);

  // create verbosity object
  mesh_ = S->GetMesh(domain_);
  dim = mesh_->getSpaceDimension();

  ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  nfaces_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // workflow can be affected by the list of models
  auto physical_models = Teuchos::sublist(ep_list_, "physical models and assumptions");
  assumptions_.Init(*physical_models, *mesh_);

  // verbose object
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ep_list_->sublist("verbose object");
  std::string ioname = "EnergyPH";
  if (domain_ != "domain") ioname += "-" + domain_;
  vo_ = Teuchos::rcp(new VerboseObject(ioname, vlist));
}


/* ******************************************************************
* Create the physical evaluators for energy, enthalpy, thermal
* conductivity, and any sources.
****************************************************************** */
void
EnergyPressureEnthalpy_PK::Setup()
{
  enthalpy_key_ = Keys::getKey(domain_, "enthalpy");
  state_key_ = Keys::getKey(domain_, "thermodynamic_state");

  energy_key_ = Keys::getKey(domain_, "energy");
  prev_energy_key_ = Keys::getKey(domain_, "prev_energy");

  temperature_key_ = Keys::getKey(domain_, "temperature");
  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");
  ie_rock_key_ = Keys::getKey(domain_, "internal_energy_rock");

  aperture_key_ = Keys::getKey(domain_, "aperture");

  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");

  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");
  porosity_key_ = Keys::getKey(domain_, "porosity");
  pressure_key_ = Keys::getKey(domain_, "pressure");
  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid");
  bcs_flow_key_ = Keys::getKey(domain_, "bcs_flow");
  bcs_enthalpy_key_ = Keys::getKey(domain_, "bcs_enthalpy");

  // require primary state variables
  std::vector<std::string> names({ "cell" });
  std::vector<AmanziMesh::Entity_kind> locations({ AmanziMesh::Entity_kind::CELL });
  std::vector<int> ndofs(1, 1);

  auto list1 = Teuchos::sublist(ep_list_, "operators", true);
  auto list2 = Teuchos::sublist(list1, "diffusion operator", true);
  auto list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  if (name != "fv: default" && name != "nlfv: default") {
    names.push_back("face");
    locations.push_back(AmanziMesh::Entity_kind::FACE);
    ndofs.push_back(1);
  } else {
    names.push_back("boundary_face");
    locations.push_back(AmanziMesh::Entity_kind::BOUNDARY_FACE);
    ndofs.push_back(1);
  }

  S_->Require<CV_t, CVS_t>(enthalpy_key_, Tags::DEFAULT)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponents(names, locations, ndofs);

  enthalpy_eval_ = Teuchos::rcp_static_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(enthalpy_key_, Tags::DEFAULT));

  // conserved quantity from the last timestep.
  if (!S_->HasRecord(prev_energy_key_)) {
    S_->Require<CV_t, CVS_t>(prev_energy_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->GetRecordW(prev_energy_key_, passwd_).set_io_vis(false);
  }

  // other fields
  // thermodynamics
  if (!S_->HasRecord(state_key_)) {
    S_->Require<CV_t, CVS_t>(state_key_, Tags::DEFAULT, state_key_, Evaluators::TS_names)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, Evaluators::TS_t_size);

    Teuchos::ParameterList elist(state_key_);
    elist.set<std::string>("tag", "");

    auto eval = Teuchos::rcp(new Evaluators::ThermodynamicStateEvaluator(elist));
    S_->SetEvaluator(state_key_, Tags::DEFAULT, eval);
  }

  // -- temperature
  if (!S_->HasRecord(temperature_key_)) {
    S_->Require<CV_t, CVS_t>(temperature_key_, Tags::DEFAULT, temperature_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    Teuchos::ParameterList elist(temperature_key_);
    elist.set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new Evaluators::TemperatureEvaluator(elist));
    S_->SetEvaluator(temperature_key_, Tags::DEFAULT, eval);
  }

  // -- rock internal energy
  if (!S_->HasRecord(ie_rock_key_)) {
    S_->Require<CV_t, CVS_t>(ie_rock_key_, Tags::DEFAULT, ie_rock_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
    S_->RequireEvaluator(ie_rock_key_, Tags::DEFAULT);
  }

  // -- densities
  if (!S_->HasRecord(mol_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    Teuchos::ParameterList elist(mol_density_liquid_key_);
    elist.set<std::string>("tag", "");

    auto eval = Teuchos::rcp(new Evaluators::DensityEvaluator(elist));
    S_->SetEvaluator(mol_density_liquid_key_, Tags::DEFAULT, eval);

    if (S_->GetEvaluator(mol_density_liquid_key_)
        .IsDifferentiableWRT(*S_, enthalpy_key_, Tags::DEFAULT)) {
      S_->RequireDerivative<CV_t, CVS_t>(mol_density_liquid_key_,
                                         Tags::DEFAULT,
                                         enthalpy_key_,
                                         Tags::DEFAULT,
                                         mol_density_liquid_key_).SetGhosted();
    }
  }

  if (!S_->HasRecord(mass_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT, mass_density_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    Teuchos::ParameterList elist(mass_density_liquid_key_);
    elist.set<std::string>("tag", "");

    auto eval = Teuchos::rcp(new Evaluators::DensityEvaluator(elist));
    S_->SetEvaluator(mass_density_liquid_key_, Tags::DEFAULT, eval);

    if (S_->GetEvaluator(mass_density_liquid_key_)
        .IsDifferentiableWRT(*S_, enthalpy_key_, Tags::DEFAULT)) {
      S_->RequireDerivative<CV_t, CVS_t>(mass_density_liquid_key_,
                                         Tags::DEFAULT,
                                         enthalpy_key_,
                                         Tags::DEFAULT,
                                         mass_density_liquid_key_).SetGhosted();
    }
  }

  // -- molar flow rates as a regular field
  if (!S_->HasRecord(mol_flowrate_key_)) {
    CompositeVectorSpace cvs;
    if (assumptions_.flow_on_manifold) {
      cvs = *Operators::CreateManifoldCVS(mesh_);
    } else {
      cvs.SetMesh(mesh_)->SetGhosted(true)->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    }

    *S_->Require<CV_t, CVS_t>(mol_flowrate_key_, Tags::DEFAULT, passwd_)
       .SetMesh(mesh_)
       ->SetGhosted(true) = cvs;
    AddDefaultPrimaryEvaluator(S_, mol_flowrate_key_, Tags::DEFAULT);
  }

  // if flow is missing, we need typical flow fields
  // -- pressure
  if (!S_->HasRecord(pressure_key_)) {
    auto cvs = S_->Require<CV_t, CVS_t>(enthalpy_key_, Tags::DEFAULT);
    *S_->Require<CV_t, CVS_t>(pressure_key_, Tags::DEFAULT, pressure_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true) = cvs;
    AddDefaultPrimaryEvaluator(S_, pressure_key_);
  }

  // -- porosity
  if (!S_->HasRecord(porosity_key_)) {
    S_->Require<CV_t, CVS_t>(porosity_key_, Tags::DEFAULT, porosity_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(porosity_key_, Tags::DEFAULT);
  }

  double molar_mass = S_->ICList().sublist("const_fluid_molar_mass").get<double>("value");

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  if (!S_->HasRecord(energy_key_)) {
    S_->Require<CV_t, CVS_t>(energy_key_, Tags::DEFAULT, energy_key_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList elist(energy_key_);
    elist.set<std::string>("tag", "");
    auto ee = Teuchos::rcp(new TotalEnergyEvaluatorPH(elist));
    S_->SetEvaluator(energy_key_, Tags::DEFAULT, ee);

    S_->RequireDerivative<CV_t, CVS_t>(
        energy_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT, energy_key_)
      .SetGhosted();

    S_->RequireDerivative<CV_t, CVS_t>(
        energy_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, energy_key_)
      .SetGhosted();
  }

  // -- thermal conductivity
  if (!S_->HasRecord(conductivity_key_)) {
    S_->Require<CV_t, CVS_t>(conductivity_key_, Tags::DEFAULT, conductivity_key_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList elist(conductivity_key_);
    elist.set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new Evaluators::ThermalConductivityEvaluator(elist));
    S_->SetEvaluator(conductivity_key_, Tags::DEFAULT, eval);
  }

  // -- viscosity
  if (!S_->HasRecord(viscosity_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(viscosity_liquid_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    Teuchos::ParameterList elist(viscosity_liquid_key_);
    elist.set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new Evaluators::ViscosityEvaluator(elist));
    S_->SetEvaluator(viscosity_liquid_key_, Tags::DEFAULT, eval);

    S_->RequireDerivative<CV_t, CVS_t>(
        viscosity_liquid_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetGhosted();

    S_->RequireDerivative<CV_t, CVS_t>(
        viscosity_liquid_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetGhosted();
  }

  // boundary conditions
  S_->Require<Operators::BCs, Operators::BCs>(bcs_enthalpy_key_, Tags::DEFAULT, "state")
    .SetMesh(mesh_)
    ->SetKind(AmanziMesh::Entity_kind::FACE)
    ->SetType(WhetStone::DOF_Type::SCALAR);
  S_->GetRecordW(bcs_enthalpy_key_, "state").set_initialized();

  // set units
  S_->GetRecordSetW(temperature_key_).set_units("K");
  S_->GetRecordSetW(mol_flowrate_key_).set_units("mol/s");
  S_->GetRecordSetW(energy_key_).set_units("J/m^3");
  S_->GetRecordSetW(enthalpy_key_).set_units("J/mol");

  if (assumptions_.flow_on_manifold) {
    S_->GetRecordSetW(aperture_key_).set_units("m");
  }
}


/* ******************************************************************
* Basic initialization of energy classes.
****************************************************************** */
void
EnergyPressureEnthalpy_PK::Initialize()
{
  // Call the base class initialize.
  Energy_PK::Initialize();

  // Create pointers to the primary field, temperature
  solution = S_->GetPtrW<CV_t>(enthalpy_key_, Tags::DEFAULT, passwd_);
  soln_->SetData(solution);

  // Create local evaluators. Initialize local fields.
  enthalpy_eval_->SetChanged();
  InitializeFields_();
  InitializeEnthalpy_();

  // initialize operators: diffusion and advection
  Teuchos::ParameterList tmp_list = ep_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  Operators::PDE_DiffusionFactory opfactory;
  Operators::PDE_AdvectionUpwindFactory opfactory_adv;

  auto op_bc_enth = S_->GetPtrW<Operators::BCs>(bcs_enthalpy_key_, Tags::DEFAULT, "state");
  op_matrix_diff_ = opfactory.Create(oplist_matrix, mesh_, op_bc_enth);
  op_matrix_diff_->SetBCs(op_bc_enth, op_bc_enth);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_matrix_->Init();

  Teuchos::RCP<Operators::BCs> bc_trial;
  if (S_->HasRecord(bcs_flow_key_)) {
    bc_trial = S_->GetPtrW<Operators::BCs>(bcs_flow_key_, Tags::DEFAULT, "state");
  } else {
    bc_trial = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
    bc_trial->bc_model();  // allocate memory 
    bc_trial->bc_value();
  } 
  op_matrix_diff_pres_ = opfactory.Create(oplist_matrix, mesh_, op_bc_enth);
  op_matrix_diff_pres_->global_operator()->Init();
  op_matrix_diff_pres_->SetBCs(bc_trial, op_bc_enth);

  Teuchos::ParameterList oplist_adv = ep_list_->sublist("operators").sublist("advection operator");
  op_matrix_advection_ = opfactory_adv.Create(oplist_adv, mesh_);

  const CompositeVector& flux = *S_->GetPtr<CV_t>(mol_flowrate_key_, Tags::DEFAULT);
  op_matrix_advection_->Setup(flux);
  op_matrix_advection_->SetBCs(op_bc_enth, op_bc_enth);
  op_advection_ = op_matrix_advection_->global_operator();

  // initialize operators: diffusion + advection + accumulation + newton_correction
  op_preconditioner_diff_ = opfactory.Create(oplist_pc, mesh_, op_bc_enth);
  op_preconditioner_diff_->SetBCs(op_bc_enth, op_bc_enth);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_preconditioner_->Init();

  op_preconditioner_advection_ = opfactory_adv.Create(oplist_adv, op_preconditioner_);
  op_preconditioner_advection_->SetBCs(op_bc_enth, op_bc_enth);

  op_acc_ = Teuchos::rcp(
    new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op_preconditioner_));

  op_preconditioner_adv_enth_ = opfactory_adv.Create(oplist_adv, op_preconditioner_);
  op_preconditioner_adv_enth_->SetBCs(op_bc_enth, op_bc_enth);

  // initialize preconditioner
  op_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  op_preconditioner_->set_inverse_parameters(pc_name, *preconditioner_list_);
  op_preconditioner_->InitializeInverse();

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  if (ti_method_name == "BDF1") {
    Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

    if (!bdf1_list.isSublist("verbose object"))
      bdf1_list.sublist("verbose object") = ep_list_->sublist("verbose object");

    bdf1_dae_ = Teuchos::rcp(
      new BDF1_TI<TreeVector, TreeVectorSpace>("BDF1", bdf1_list, *this, soln_->get_map(), S_));
  }

  // initialize boundary conditions
  double t_ini = S_->get_time();
  auto& temperature = S_->GetW<CV_t>(temperature_key_, temperature_key_);
  UpdateSourceBoundaryData(t_ini, t_ini, temperature);

  // output of initialization summary
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "temperature BC assigned to " << dirichlet_bc_faces_ << " faces\n"
               << "default (zero-gradient) BC assigned to " << missed_bc_faces_ << " faces\n\n" 
               << "solution vector: ";
    solution->Print(*vo_->os(), false);
    *vo_->os() << "matrix: " << my_operator(Operators::OPERATOR_MATRIX)->PrintDiagnostics()
               << std::endl
               << "preconditioner: "
               << my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->PrintDiagnostics()
               << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl;
  }

  if (dirichlet_bc_faces_ == 0 && domain_ == "domain" && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }
}


/* *******************************************************************
* Initialize additional fields
******************************************************************* */
void 
EnergyPressureEnthalpy_PK::InitializeEnthalpy_()
{
  S_->GetEvaluator(temperature_key_).Update(*S_, "energy");

  const auto& T_c = *S_->Get<CV_t>(temperature_key_).ViewComponent("cell");
  const auto& p_c = *S_->Get<CV_t>(pressure_key_).ViewComponent("cell");
  auto& h_c = *S_->Get<CV_t>(enthalpy_key_).ViewComponent("cell");

  Teuchos::ParameterList plist;
  auto eos = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));

  int count = h_c.MyLength();
  for (int i = 0; i != count; ++i) {
    double pMPa = p_c[0][i] * 1.0e-6;
    h_c[0][i] = cfactor * (eos->ThermodynamicsPT(pMPa, T_c[0][i])).h;
  }

  DeriveFaceValuesFromCellValues(S_->GetW<CV_t>(enthalpy_key_, Tags::DEFAULT, passwd_));
}


/* *******************************************************************
* Performs one timestep of size t_new - t_old
******************************************************************* */
bool
EnergyPressureEnthalpy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // save a copy of primary and conservative fields
  std::vector<std::string> fields({ enthalpy_key_, energy_key_ });
  archive_ = Teuchos::rcp(new StateArchive(S_, vo_));
  archive_->Add(fields, Tags::DEFAULT);

  // initialization
  if (num_itrs_ == 0) {
    Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae_->SetInitialState(t_old, soln_, udot);

    UpdatePreconditioner(t_old, soln_, dt_);
    num_itrs_++;
  }

  // trying to make a step
  bool failed = bdf1_dae_->AdvanceStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;
    archive_->Restore("");
    enthalpy_eval_->SetChanged();
  }

  dt_ = dt_next_;
  return failed;
}


/* ******************************************************************
* TBW
****************************************************************** */
void
EnergyPressureEnthalpy_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  // commit solution to time history
  if (bdf1_dae_.get()) bdf1_dae_->CommitSolution(t_new - t_old, soln_);
  enthalpy_eval_->SetChanged();

  num_itrs_++;

  // update previous fields
  std::vector<std::string> fields({ energy_key_ });
  StateArchive archive_tmp(S_, vo_);
  archive_tmp.CopyFieldsToPrevFields(fields, "", false);
}


/* ******************************************************************
* Restore state to the previous step
****************************************************************** */
void
EnergyPressureEnthalpy_PK::FailStep(double t_old, double t_new, const Tag& tag)
{
  archive_->Restore("");
  enthalpy_eval_->SetChanged();
}


/* ******************************************************************
* New implementation of additional boundary conditions.
****************************************************************** */
void EnergyPressureEnthalpy_PK::ComputeSecondaryBCs() 
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  Teuchos::RCP<const Epetra_MultiVector> p_f, p_bf;
  Teuchos::RCP<Epetra_MultiVector> h_f, h_bf;

  const auto& p = S_->Get<CV_t>(pressure_key_);
  auto& h = S_->GetW<CV_t>(enthalpy_key_, Tags::DEFAULT, passwd_);

  if (p.HasComponent("face")) p_f = p.ViewComponent("face", true);
  else p_bf = p.ViewComponent("boundary_face", true);

  if (h.HasComponent("face")) h_f = h.ViewComponent("face", true);
  else h_bf = h.ViewComponent("boundary_face", true);

  auto& T_bf = *S_->GetW<CV_t>(temperature_key_, Tags::DEFAULT, temperature_key_).ViewComponent("boundary_face", true);

  Teuchos::ParameterList plist;
  auto eos = Teuchos::rcp(new AmanziEOS::IAPWS97(plist));

  // populate BCs
  auto op_bc_enth = S_->GetPtrW<Operators::BCs>(bcs_enthalpy_key_, Tags::DEFAULT, "state");
  std::vector<int>& bc_model_enth = op_bc_enth->bc_model();
  std::vector<double>& bc_value_enth = op_bc_enth->bc_value();

  for (int n = 0; n < bc_model.size(); ++n) {
    bc_model_enth[n] = Operators::OPERATOR_BC_NONE;
    bc_value_enth[n] = 0.0;
  }

  int nbfaces = mesh_->getNumEntities(AmanziMesh::Entity_kind::BOUNDARY_FACE, AmanziMesh::Parallel_kind::OWNED);
  double pMPa;
  for (int bf = 0; bf < nbfaces; ++bf) {
    int f = getBoundaryFaceFace(*mesh_, bf);
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      if (p_f.get()) pMPa = (*p_f)[0][f] * 1e-6;
      else pMPa = (*p_bf)[0][bf] * 1e-6;

      bc_model_enth[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value_enth[f] = cfactor * eos->ThermodynamicsPT(pMPa, bc_value[f]).h;

      if (p_f.get()) (*h_f)[0][f] = bc_value_enth[f];
      else (*h_bf)[0][bf] = bc_value_enth[f];

      T_bf[0][bf] = bc_value[f];
    }
    else if (bc_model[f] == Operators::OPERATOR_BC_TOTAL_FLUX) {
      bc_model_enth[f] = Operators::OPERATOR_BC_TOTAL_FLUX;
      bc_value_enth[f] = bc_value[f];
    }
  }
}

} // namespace Energy
} // namespace Amanzi
