/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the Energy component of the Amanzi code.

  Process kernel for the single-phase energy equation.
*/

// Amanzi
#include "PDE_DiffusionFactory.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "UpwindFactory.hh"

// Amanzi::Energy
#include "EnergyOnePhase_PK.hh"
#include "EnthalpyEvaluator.hh"
#include "IEMEvaluator.hh"
#include "StateArchive.hh"
#include "TCMEvaluator_OnePhase.hh"
#include "TotalEnergyEvaluator.hh"

namespace Amanzi {
namespace Energy {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Default constructor.
****************************************************************** */
EnergyOnePhase_PK::EnergyOnePhase_PK(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln),
    Energy_PK(pk_tree, glist, S, soln),
    soln_(soln),
    num_itrs_(0),
    dt_(0.0)
{
  // verbose object
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ep_list_->sublist("verbose object");
  std::string ioname = "Energy1Phase";
  if (domain_ != "domain") ioname += "-" + domain_;
  vo_ = Teuchos::rcp(new VerboseObject(ioname, vlist));
}


/* ******************************************************************
* Create the physical evaluators for energy, enthalpy, thermal
* conductivity, and any sources.
****************************************************************** */
void
EnergyOnePhase_PK::Setup()
{
  // basic class setup
  Energy_PK::Setup();

  double molar_mass = S_->ICList().sublist("const_fluid_molar_mass").get<double>("value");

  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  if (!S_->HasRecord(energy_key_)) {
    S_->Require<CV_t, CVS_t>(energy_key_, Tags::DEFAULT, energy_key_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList elist = ep_list_->sublist("energy evaluator");
    elist.set<std::string>("energy key", energy_key_)
      .set<std::string>("particle density key", particle_density_key_)
      .set<std::string>("internal energy rock key", ie_rock_key_)
      .set<bool>("vapor diffusion", false)
      .set<double>("liquid molar mass", molar_mass)
      .set<std::string>("tag", "");
    if (flow_on_manifold_) elist.set<std::string>("aperture key", aperture_key_);

    elist.setName(energy_key_);
    auto ee = Teuchos::rcp(new TotalEnergyEvaluator(elist));
    S_->SetEvaluator(energy_key_, Tags::DEFAULT, ee);

    S_->RequireDerivative<CV_t, CVS_t>(
        energy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, energy_key_)
      .SetGhosted();
  }

  // -- advection of enthalpy
  if (!S_->HasRecord(enthalpy_key_)) {
    S_->Require<CV_t, CVS_t>(enthalpy_key_, Tags::DEFAULT, enthalpy_key_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    Teuchos::ParameterList elist = ep_list_->sublist("enthalpy evaluator");
    elist.set("enthalpy key", enthalpy_key_)
      .set<double>("liquid molar mass", molar_mass)
      .set<std::string>("tag", "");
    elist.setName(enthalpy_key_);
    auto enth = Teuchos::rcp(new EnthalpyEvaluator(elist));
    S_->SetEvaluator(enthalpy_key_, Tags::DEFAULT, enth);

    S_->RequireDerivative<CV_t, CVS_t>(
        enthalpy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, enthalpy_key_)
      .SetGhosted();
  }

  // -- thermal conductivity
  if (!S_->HasRecord(conductivity_key_)) {
    S_->Require<CV_t, CVS_t>(conductivity_key_, Tags::DEFAULT, conductivity_key_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList elist = ep_list_->sublist("thermal conductivity evaluator");
    elist.set("thermal conductivity key", conductivity_key_).set<std::string>("tag", "");
    elist.setName(conductivity_key_);

    auto tcm = Teuchos::rcp(new TCMEvaluator_OnePhase(mesh_, elist));
    S_->SetEvaluator(conductivity_key_, Tags::DEFAULT, tcm);
  }

  // set units
  S_->GetRecordSetW(energy_key_).set_units("J/m^3");
  S_->GetRecordSetW(enthalpy_key_).set_units("J/mol");
}


/* ******************************************************************
* Initialize the needed models to plug in enthalpy.
****************************************************************** */
void
EnergyOnePhase_PK::Initialize()
{
  // Call the base class initialize.
  Energy_PK::Initialize();

  // Create pointers to the primary field, temperature
  solution = S_->GetPtrW<CV_t>(temperature_key_, Tags::DEFAULT, passwd_);
  soln_->SetData(solution);

  // Create local evaluators. Initialize local fields.
  InitializeFields_();

  // initialize independent operators: diffusion and advection
  Teuchos::ParameterList tmp_list = ep_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  Operators::PDE_DiffusionFactory opfactory;
  Operators::PDE_AdvectionUpwindFactory opfactory_adv;

  op_matrix_diff_ = opfactory.Create(oplist_matrix, mesh_, op_bc_);
  op_matrix_diff_->SetBCs(op_bc_, op_bc_);
  op_matrix_ = op_matrix_diff_->global_operator();
  op_matrix_->Init();

  Teuchos::ParameterList oplist_adv = ep_list_->sublist("operators").sublist("advection operator");
  op_matrix_advection_ = opfactory_adv.Create(oplist_adv, mesh_);

  const CompositeVector& flux = *S_->GetPtr<CV_t>(mol_flowrate_key_, Tags::DEFAULT);
  op_matrix_advection_->Setup(flux);
  op_matrix_advection_->SetBCs(op_bc_enth_, op_bc_enth_);
  op_advection_ = op_matrix_advection_->global_operator();

  // initialize coupled operators: diffusion + advection + accumulation
  op_preconditioner_diff_ = opfactory.Create(oplist_pc, mesh_, op_bc_);
  op_preconditioner_diff_->SetBCs(op_bc_, op_bc_);
  op_preconditioner_ = op_preconditioner_diff_->global_operator();
  op_preconditioner_->Init();

  // optional upwinding of conductivity
  if (tmp_list.isSublist("conductivity")) {
    Operators::UpwindFactory upw_factory;
    upwind_ = upw_factory.Create(mesh_, tmp_list.sublist("conductivity"));
    upw_conductivity_ = Teuchos::rcp(new CompositeVector(*upwind_->Map()));

    op_matrix_diff_->SetScalarCoefficient(upw_conductivity_, Teuchos::null);
    op_preconditioner_diff_->SetScalarCoefficient(upw_conductivity_, Teuchos::null);
  } else {
    op_matrix_diff_->SetScalarCoefficient(S_->GetPtr<CV_t>(conductivity_gen_key_, Tags::DEFAULT),
                                          Teuchos::null);
    op_preconditioner_diff_->SetScalarCoefficient(
      S_->GetPtr<CV_t>(conductivity_gen_key_, Tags::DEFAULT), Teuchos::null);
  }

  op_acc_ = Teuchos::rcp(
    new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op_preconditioner_));
  if (prec_include_enthalpy_) {
    op_preconditioner_advection_ = opfactory_adv.Create(oplist_adv, op_preconditioner_);
    op_preconditioner_advection_->SetBCs(op_bc_enth_, op_bc_enth_);
  }

  // initialize preconditioner
  op_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  std::string name = ti_list_->get<std::string>("preconditioner");
  Teuchos::ParameterList slist = preconditioner_list_->sublist(name);
  op_preconditioner_->set_inverse_parameters(slist);
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
  auto& temperature = S_->GetW<CV_t>(temperature_key_, passwd_);
  UpdateSourceBoundaryData(t_ini, t_ini, temperature);

  // output of initialization summary
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "temperature BC assigned to " << dirichlet_bc_faces_ << " faces\n\n"
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


/* ****************************************************************
* This completes initialization of missed fields in the state.
**************************************************************** */
void
EnergyOnePhase_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (S_->HasRecord(prev_energy_key_)) {
    if (!S_->GetRecord(prev_energy_key_).initialized()) {
      temperature_eval_->SetChanged();
      S_->GetEvaluator(energy_key_).Update(*S_, passwd_);

      const auto& e1 = S_->Get<CV_t>(energy_key_);
      auto& e0 = S_->GetW<CV_t>(prev_energy_key_, passwd_);
      e0 = e1;

      S_->GetRecordW(prev_energy_key_, passwd_).set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
        *vo_->os() << "initialized prev_energy to previous energy" << std::endl;
    }
  }
}


/* *******************************************************************
* Performs one timestep of size dt_ either for steady-state or
* transient sumulation.
******************************************************************* */
bool
EnergyOnePhase_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  // save a copy of primary and conservative fields
  std::vector<std::string> fields({ temperature_key_, energy_key_ });
  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

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

    archive.Restore("");
    temperature_eval_->SetChanged();
    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae_->CommitSolution(dt_, soln_);
  temperature_eval_->SetChanged();

  num_itrs_++;
  dt_ = dt_next_;

  return failed;
}


/* ******************************************************************
* TBW
****************************************************************** */
void
EnergyOnePhase_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  dt_ = dt_next_;

  // update previous fields
  std::vector<std::string> fields({ energy_key_ });
  StateArchive archive(S_, vo_);
  archive.CopyFieldsToPrevFields(fields, "", false);
}

} // namespace Energy
} // namespace Amanzi
