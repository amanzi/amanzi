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

#include "InverseFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PorosityEvaluator.hh"
#include "PDE_ElasticityFactory.hh"
#include "StateArchive.hh"
#include "PK_Helpers.hh"

#include "HydrostaticStressEvaluator.hh"
#include "MechanicsElasticity_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
MechanicsElasticity_PK::MechanicsElasticity_PK(Teuchos::ParameterList& pk_tree,
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
  requireAtNext(displacement_key_, Tags::DEFAULT, *S_, passwd_);

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
MechanicsElasticity_PK::Setup()
{
  Mechanics_PK::Setup();
}


/* ******************************************************************
* Function goes through flow parameter list and initializes various
* objects including those created during the setup step.
****************************************************************** */
void
MechanicsElasticity_PK::Initialize()
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
  Operators::PDE_ElasticityFactory opfactory;

  // -- create elastic block
  auto tmp1 = ec_list_->sublist("operators").sublist("elasticity operator");
  op_matrix_elas_ = opfactory.Create(tmp1, mesh_);

  op_matrix_ = op_matrix_elas_->global_operator();

  // -- extensions: The undrained split method add anotehr operator which has
  //    the grad-div structure. It is critical that it uses a separate global
  //    operator pointer. Its local matrices are shared with the original
  //    physics operator.
  if (split_undrained_) {
    std::string method = tmp1.sublist("schema").get<std::string>("method");
    tmp1.sublist("schema").set<std::string>("method", method + " graddiv");
    op_matrix_graddiv_ = opfactory.Create(tmp1, mesh_);
    op_matrix_->OpPushBack(op_matrix_graddiv_->local_op());
  }

  // Create BC objects
  InitializeBCs();

  // Populate matrix and preconditioner
  // -- setup phase
  const auto& E = S_->GetPtr<CV_t>(young_modulus_key_, Tags::DEFAULT);
  const auto& nu = S_->GetPtr<CV_t>(poisson_ratio_key_, Tags::DEFAULT);

  op_matrix_->Init();
  op_matrix_elas_->SetTensorCoefficientEnu(E, nu);

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
MechanicsElasticity_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  StateArchive archive(S_, vo_);
  archive.Add({ displacement_key_ }, Tags::DEFAULT);

  UpdateSourceBoundaryData(t_old, t_new);

  auto op = op_matrix_elas_->global_operator();
  op->Init();

  // add external forces
  auto rhs = op_matrix_->rhs();
  if (use_gravity_) AddGravityTerm(*rhs);
  if (poroelasticity_) AddPressureGradient(*rhs);
  if (thermoelasticity_) AddTemperatureGradient(*rhs);

  // update the matrix = preconditioner
  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  // extensions
  if (split_undrained_) {
    auto u = S_->Get<CV_t>(displacement_key_, Tags::DEFAULT);
    op_matrix_graddiv_->SetScalarCoefficient(
      S_->Get<CV_t>(undrained_split_coef_key_, Tags::DEFAULT));
    op_matrix_graddiv_->UpdateMatrices();
    op_matrix_graddiv_->ApplyBCs(false, true, false);
    op_matrix_graddiv_->global_operator()->Apply(u, *rhs, 1.0);
  }

  // solver the problem
  op->ComputeInverse();
  int ierr = op->ApplyInverse(*rhs, *solution_);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "elasticity solver (PCG): ||r||_H=" << op->residual() << " itr=" << op->num_itrs()
               << " code=" << op->returned_code() << std::endl;
  }

  if (ierr != 0) {
    archive.Restore("");
    dt_ /= 2;
    return true;
  }

  dt_ *= 2;
  eval_->SetChanged();
  return false;
}


/* *******************************************************************
* Performs one timestep from time t_old to time t_new either for
* steady-state or transient simulation.
******************************************************************* */
void
MechanicsElasticity_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetEvaluator(hydrostatic_stress_key_).Update(*S_, "mechanics");
  S_->GetEvaluator(vol_strain_key_).Update(*S_, "mechanics");
}


/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
Teuchos::RCP<Operators::Operator>
MechanicsElasticity_PK::my_operator(const Operators::OperatorType& type)
{
  return op_matrix_;
}

} // namespace Mechanics
} // namespace Amanzi
