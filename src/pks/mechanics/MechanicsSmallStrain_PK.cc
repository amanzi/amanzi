/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

*/

#include <vector>

// Amanzi
#include "InverseFactory.hh"
#include "PDE_ElasticityFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PorosityEvaluator.hh"
#include "StateArchive.hh"
#include "StateHelpers.hh"

#include "HydrostaticStressEvaluator.hh"
#include "MechanicsSmallStrain_PK.hh"
#include "ShearStrainEvaluator.hh"
#include "SSMEvaluator.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
MechanicsSmallStrain_PK::MechanicsSmallStrain_PK(Teuchos::ParameterList& pk_tree,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                                 const Teuchos::RCP<State>& S,
                                                 const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln), Mechanics_PK(pk_tree, glist, S, soln), soln_(soln)
{
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ec_list_ = Teuchos::sublist(pk_list, name_, true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ec_list_, "time integrator", true);

  // domain and primary evaluators
  domain_ = ec_list_->get<std::string>("domain name", "domain");
  displacement_key_ = Keys::getKey(domain_, "displacement");
  AddDefaultPrimaryEvaluator(S_, displacement_key_);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ec_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("MechanicsElasticity", vlist));
}


/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void
MechanicsSmallStrain_PK::Setup()
{
  Mechanics_PK::Setup();

  shear_strain_key_ = Keys::getKey(domain_, "shear_strain");
  shear_modulus_key_ = Keys::getKey(domain_, "shear_modulus");
  bulk_modulus_key_ = Keys::getKey(domain_, "bulk_modulus");

  if (!S_->HasRecord(shear_strain_key_)) {
    auto elist = RequireFieldForEvaluator(*S_, shear_strain_key_);

    eval_shear_strain_ = Teuchos::rcp(new ShearStrainEvaluator(elist));
    S_->SetEvaluator(shear_strain_key_, Tags::DEFAULT, eval_shear_strain_);
  }

  auto ssm_list = Teuchos::sublist(ec_list_, "small strain models", true);
  auto ssm = CreateSSMPartition(mesh_, ssm_list);

  if (!S_->HasRecord(shear_modulus_key_)) {
    auto elist = RequireFieldForEvaluator(*S_, shear_modulus_key_);

    auto eval = Teuchos::rcp(new SSMEvaluator(elist, ssm));
    S_->SetEvaluator(shear_modulus_key_, Tags::DEFAULT, eval);
  }

  S_->Require<CV_t, CVS_t>(bulk_modulus_key_, Tags::DEFAULT, bulk_modulus_key_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(bulk_modulus_key_, Tags::DEFAULT);
}


/* ******************************************************************
* Function goes through flow parameter list and initializes various
* objects including those created during the setup step.
****************************************************************** */
void
MechanicsSmallStrain_PK::Initialize()
{
  Mechanics_PK::Initialize();

  // -- miscalleneous defaults
  num_itrs_ = 0;
  double t_ini = S_->get_time();

  // -- control parameters
  auto physical_models = Teuchos::sublist(ec_list_, "physical models and assumptions");
  use_gravity_ = physical_models->get<bool>("use gravity");
  thermoelasticity_ = physical_models->get<bool>("thermoelasticity", false);

  split_undrained_ = physical_models->get<bool>("biot scheme: undrained split", false);
  split_fixed_stress_ = physical_models->get<bool>("biot scheme: fixed stress split", false);
  poroelasticity_ = split_undrained_ || split_fixed_stress_;

  // Create pointers to the primary flow field displacement.
  solution_ = S_->GetPtrW<CV_t>(displacement_key_, Tags::DEFAULT, passwd_);
  soln_->SetData(solution_);

  // Initialize time integrator.
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  if (ti_method_name == "BDF1") {
    Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");
    bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>("BDF1", bdf1_list, *this, soln_->get_map(), S_));
  }

  // Initialize matrix and preconditioner
  // -- create elastic block
  auto tmp1 = ec_list_->sublist("operators").sublist("elasticity operator");
  Operators::PDE_ElasticityFactory factory;
  op_matrix_elas_ = factory.Create(tmp1, mesh_);

  op_matrix_ = op_matrix_elas_->global_operator();
  // For a quasi-static problem, we update PC every non-linear iteration.
  // op_preconditioner_elas_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));
  // op_preconditioner_elas_->Init(tmp1);

  // -- extensions: The undrained split method add another operator which has
  //    the grad-div structure. It is critical that it uses a separate global
  //    operator pointer. Its local matrices are shared with the original
  //    physics operator.
  if (split_undrained_) {
    std::string method = tmp1.sublist("schema").get<std::string>("method");
    tmp1.sublist("schema").set<std::string>("method", method + " graddiv");
    op_matrix_graddiv_ = factory.Create(tmp1, mesh_);

    op_matrix_->OpPushBack(op_matrix_graddiv_->local_op());
  }

  // Create BC objects
  InitializeBCs();

  // Populate matrix and preconditioner
  // -- setup phase
  const auto& G = S_->GetPtr<CV_t>(shear_modulus_key_, Tags::DEFAULT);
  const auto& K = S_->GetPtr<CV_t>(bulk_modulus_key_, Tags::DEFAULT);

  op_matrix_->Init();
  op_matrix_elas_->SetTensorCoefficientGK(G, K);

  for (auto bc : op_bcs_) op_matrix_elas_->AddBCs(bc, bc);

  UpdateSourceBoundaryData(t_ini, t_ini);
  op_matrix_elas_->EnforceBCs(*solution_);

  // -- assemble phase
  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  // -- preconditioned linear solver for alternative solution strategy
  AMANZI_ASSERT(ti_list_->isParameter("linear solver"));
  std::string solver_name = ti_list_->get<std::string>("linear solver");

  AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  op_matrix_->set_inverse_parameters(
    pc_name, *preconditioner_list_, solver_name, *linear_solver_list_, true);

  op_matrix_->InitializeInverse();

  // we set up operators and can trigger re-initialization of stress
  eval_hydro_stress_->set_op(op_matrix_elas_);
  eval_vol_strain_->set_op(op_matrix_elas_);
  eval_shear_strain_->set_op(op_matrix_elas_);
  eval_->SetChanged();

  // summary of initialization
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << " TI:\"" << ti_method_name.c_str() << "\"" << std::endl
               << "matrix: " << op_matrix_->PrintDiagnostics() << std::endl;

    *vo_->os() << "displacement BC assigned to " << dirichlet_bc_ << " nodes" << std::endl;

    *vo_->os() << vo_->color("green") << std::endl
               << "Initialization of PK is complete, T=" << units_.OutputTime(S_->get_time())
               << vo_->reset() << std::endl
               << std::endl;
  }
}


/* *******************************************************************
* Performs one timestep from time t_old to time t_new either for
* steady-state or transient simulation.
******************************************************************* */
bool
MechanicsSmallStrain_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  StateArchive archive(S_, vo_);
  archive.Add({ displacement_key_ }, Tags::DEFAULT);

  // initialization
  if (num_itrs_ == 0) {
    auto udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae_->SetInitialState(t_old, soln_, udot);

    UpdatePreconditioner(t_old, soln_, dt_);
    num_itrs_++;
  }

  // trying to make a step
  bool failed(false);
  failed = bdf1_dae_->AdvanceStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;

    archive.Restore("");
    eval_->SetChanged();

    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae_->CommitSolution(dt_, soln_);
  eval_->SetChanged();
  num_itrs_++;

  return false;
}


/* *******************************************************************
* Updates after successful timesteps
******************************************************************* */
void
MechanicsSmallStrain_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetEvaluator(hydrostatic_stress_key_).Update(*S_, "mechanics");
  S_->GetEvaluator(vol_strain_key_).Update(*S_, "mechanics");
}

} // namespace Mechanics
} // namespace Amanzi
