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
  : soln_(soln), passwd_("")
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
MechanicsElasticity_PK::Setup()
{
  dt_ = 1e+98;
  mesh_ = S_->GetMesh();
  dim_ = mesh_->getSpaceDimension();

  displacement_key_ = Keys::getKey(domain_, "displacement");
  hydrostatic_stress_key_ = Keys::getKey(domain_, "hydrostatic_stress");
  vol_strain_key_ = Keys::getKey(domain_, "volumetric_strain");

  young_modulus_key_ = Keys::getKey(domain_, "young_modulus");
  poisson_ratio_key_ = Keys::getKey(domain_, "poisson_ratio");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");
  undrained_split_coef_key_ = Keys::getKey(domain_, "undrained_split_coef");

  // constant fields
  S_->Require<AmanziGeometry::Point>("gravity", Tags::DEFAULT, "state");

  // primary fields
  // -- displacement
  const auto& list = ec_list_->sublist("operators").sublist("elasticity operator");
  auto schema = Operators::schemaFromPList(list, mesh_);

  if (!S_->HasRecord(displacement_key_)) {
    auto cvs = Operators::cvsFromSchema(schema, mesh_, true);
    *S_->Require<CV_t, CVS_t>(displacement_key_, Tags::DEFAULT, passwd_)
       .SetMesh(mesh_)
       ->SetGhosted(true) = cvs;

    Teuchos::ParameterList elist(displacement_key_);
    elist.set<std::string>("evaluator name", displacement_key_);
    eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(displacement_key_, Tags::DEFAULT, eval_);
  }

  if (!S_->HasRecord(hydrostatic_stress_key_)) {
    S_->Require<CV_t, CVS_t>(hydrostatic_stress_key_, Tags::DEFAULT, hydrostatic_stress_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }
  {
    Teuchos::ParameterList elist(hydrostatic_stress_key_);
    elist.set<std::string>("tag", "");
    eval_hydro_stress_ = Teuchos::rcp(new HydrostaticStressEvaluator(elist));
    S_->SetEvaluator(hydrostatic_stress_key_, Tags::DEFAULT, eval_hydro_stress_);
  }

  if (!S_->HasRecord(vol_strain_key_)) {
    S_->Require<CV_t, CVS_t>(vol_strain_key_, Tags::DEFAULT, vol_strain_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1); // copy of states needs it
  }
  {
    Teuchos::ParameterList elist(vol_strain_key_);
    elist.set<std::string>("tag", "");
    // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    eval_vol_strain_ = Teuchos::rcp(new VolumetricStrainEvaluator(elist));
    S_->SetEvaluator(vol_strain_key_, Tags::DEFAULT, eval_vol_strain_);
  }

  // -- rock properties
  S_->Require<CV_t, CVS_t>(young_modulus_key_, Tags::DEFAULT, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  S_->Require<CV_t, CVS_t>(poisson_ratio_key_, Tags::DEFAULT, passwd_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  if (!S_->HasRecord(particle_density_key_)) {
    S_->Require<CV_t, CVS_t>(particle_density_key_, Tags::DEFAULT, particle_density_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(particle_density_key_, Tags::DEFAULT);
  }

  // set units
  S_->GetRecordSetW(displacement_key_).set_units("m");
  S_->GetRecordSetW(poisson_ratio_key_).set_units("-");
  S_->GetRecordSetW(young_modulus_key_).set_units("Pa");
}


/* ******************************************************************
* Function goes through flow parameter list and initializes various
* objects including those created during the setup step.
****************************************************************** */
void
MechanicsElasticity_PK::Initialize()
{
  // Initialize miscalleneous defaults.
  num_itrs_ = 0;

  // -- times
  double t_ini = S_->get_time();

  // -- mesh dimensions
  ncells_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  nfaces_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  nnodes_owned_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_wghost_ =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  // -- control parameters
  auto physical_models = Teuchos::sublist(ec_list_, "physical models and assumptions");
  use_gravity_ = physical_models->get<bool>("use gravity");
  thermoelasticity_ = physical_models->get<bool>("thermoelasticity", false);

  split_undrained_ = physical_models->get<bool>("biot scheme: undrained split", false);
  split_fixed_stress_ = physical_models->get<bool>("biot scheme: fixed stress split", false);
  poroelasticity_ = split_undrained_ || split_fixed_stress_;

  // Create verbosity object to print out initialiation statistics.
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nPK initialization started...\n";
  }

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

  // -- extensions: The undrained split method add anotehr operator which has
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
  const auto E = S_->GetPtr<CV_t>(young_modulus_key_, Tags::DEFAULT);
  const auto nu = S_->GetPtr<CV_t>(poisson_ratio_key_, Tags::DEFAULT);

  op_matrix_->Init();
  op_matrix_elas_->SetTensorCoefficient(E, nu);

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

  UpdateSourceBoundaryData_(t_ini, t_ini);
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
MechanicsElasticity_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  dt_ = t_new - t_old;

  StateArchive archive(S_, vo_);
  archive.Add({ displacement_key_ }, Tags::DEFAULT);

  UpdateSourceBoundaryData_(t_old, t_new);

  auto op = op_matrix_elas_->global_operator();
  op->Init();

  // add external forces
  auto rhs = op_matrix_->rhs();
  if (use_gravity_) AddGravityTerm_(*rhs);
  if (poroelasticity_) AddPressureGradient_(*rhs);
  if (thermoelasticity_) AddTemperatureGradient_(*rhs);

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
* Performs one time step from time t_old to time t_new either for
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
