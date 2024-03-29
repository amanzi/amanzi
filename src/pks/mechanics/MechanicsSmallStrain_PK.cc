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
#include "StateArchive.hh"
#include "StateHelpers.hh"

#include "HydrostaticStressEvaluator.hh"
#include "MechanicsSmallStrain_PK.hh"
#include "ShearModulusEvaluator.hh"

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
  : Mechanics_PK(pk_tree, glist, S, soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ec_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ec_list_, "time integrator", true);

  // domain name
  domain_ = ec_list_->get<std::string>("domain name", "domain");

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

  shear_modulus_key_ = Keys::getKey(domain_, "shear_modulus");
  bulk_modulus_key_ = Keys::getKey(domain_, "bulk_modulus");

  if (!S_->HasRecord(shear_modulus_key_)) {
    auto elist = RequireFieldForEvaluator(*S_, shear_modulus_key_);

    auto eval = Teuchos::rcp(new ShearModulusEvaluator(elist));
    S_->SetEvaluator(shear_modulus_key_, Tags::DEFAULT, eval);
  }

  S_->Require<CV_t, CVS_t>(bulk_modulus_key_, Tags::DEFAULT, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
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
    bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  }

  // Initialize matrix and preconditioner
  // -- create elastic block
  auto tmp1 = ec_list_->sublist("operators").sublist("elasticity operator");
  op_matrix_elas_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));
  op_matrix_ = op_matrix_elas_->global_operator();

  // -- extensions: The undrained split method add another operator which has
  //    the grad-div structure. It is critical that it uses a separate global
  //    operator pointer. Its local matrices are shared with the original
  //    physics operator.
  if (split_undrained_) {
    std::string method = tmp1.sublist("schema").get<std::string>("method");
    tmp1.sublist("schema").set<std::string>("method", method + " graddiv");
    op_matrix_graddiv_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));
    op_matrix_->OpPushBack(op_matrix_graddiv_->local_op());
  }

  // Create BC objects
  Teuchos::RCP<Teuchos::ParameterList> bc_list =
    Teuchos::rcp(new Teuchos::ParameterList(ec_list_->sublist("boundary conditions", true)));

  bcs_.clear();

  // -- displacement
  if (bc_list->isSublist("displacement")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("displacement");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        // nodal dofs
        auto bc =
          bc_factory.Create(spec, "no slip", AmanziMesh::NODE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("no slip");
        bc->set_type(WhetStone::DOF_Type::POINT);
        bc->set_kind(AmanziMesh::NODE);
        bcs_.push_back(bc);

        // bubble dofs
        auto bc2 =
          bc_factory.Create(spec, "no slip", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc2->set_bc_name("no slip");
        bc2->set_type(WhetStone::DOF_Type::POINT);
        bc2->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc2);
      }
    }
  }

  if (bc_list->isSublist("kinematic")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("kinematic");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        // nodal dofs
        auto bc = bc_factory.Create(
          spec, "kinematic", AmanziMesh::NODE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("kinematic");
        bc->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bc->set_kind(AmanziMesh::NODE);
        bcs_.push_back(bc);

        // bubble dofs
        auto bc2 = bc_factory.Create(
          spec, "kinematic", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc2->set_bc_name("kinematic");
        bc2->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bc2->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc2);
      }
    }
  }

  if (bc_list->isSublist("traction")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("traction");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        auto bc =
          bc_factory.Create(spec, "traction", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("traction");
        bc->set_type(WhetStone::DOF_Type::POINT);
        bc->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc);
      }
    }
  }

  if (bc_list->isSublist("normal traction")) {
    PK_DomainFunctionFactory<MechanicsBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("normal traction");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        auto bc = bc_factory.Create(
          spec, "normal traction", AmanziMesh::FACE, Teuchos::null, Tags::DEFAULT, true);
        bc->set_bc_name("normal traction");
        bc->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bc->set_kind(AmanziMesh::FACE);
        bcs_.push_back(bc);
      }
    }
  }

  // Populate matrix and preconditioner
  // -- setup phase
  const auto& G = S_->GetPtr<CV_t>(shear_modulus_key_, Tags::DEFAULT);
  const auto& K = S_->GetPtr<CV_t>(bulk_modulus_key_, Tags::DEFAULT);

  op_matrix_->Init();
  op_matrix_elas_->SetTensorCoefficientGK(G, K);

  // -- initialize boundary conditions (memory allocation)
  auto bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::POINT));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  op_bcs_.push_back(bc);

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
  Teuchos::rcp_dynamic_cast<HydrostaticStressEvaluator>(eval_hydro_stress_)
    ->set_op(op_matrix_elas_);
  Teuchos::rcp_dynamic_cast<VolumetricStrainEvaluator>(eval_vol_strain_)->set_op(op_matrix_elas_);
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
* Performs one time step from time t_old to time t_new either for
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
  failed = bdf1_dae_->TimeStep(dt_, dt_next_, soln_);
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
* Performs one time step from time t_old to time t_new either for
* steady-state or transient simulation.
******************************************************************* */
void
MechanicsSmallStrain_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  S_->GetEvaluator(hydrostatic_stress_key_).Update(*S_, "mechanics");
  S_->GetEvaluator(vol_strain_key_).Update(*S_, "mechanics");
}

} // namespace Mechanics
} // namespace Amanzi
