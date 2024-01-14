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

#include "HydrostaticStressEvaluator.hh"
#include "MechanicsElasticity_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
MechanicsElasticity_PK::MechanicsElasticity_PK(
  Teuchos::ParameterList& pk_tree,
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
  dt_desirable_ = dt_;
  dt_next_ = dt_;

  // -- mesh dimensions
  ncells_owned_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  nfaces_owned_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  nnodes_owned_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_wghost_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  // -- control parameters
  use_gravity_ = ec_list_->sublist("physical models and assumptions").get<bool>("use gravity");
  biot_model_ = ec_list_->sublist("physical models and assumptions").get<bool>("use biot model");

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
  Teuchos::ParameterList& tmp1 = ec_list_->sublist("operators").sublist("elasticity operator");
  op_matrix_elas_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));
  op_preconditioner_elas_ = Teuchos::rcp(new Operators::PDE_Elasticity(tmp1, mesh_));

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

        auto bc = bc_factory.Create(spec, "no slip", AmanziMesh::NODE, Teuchos::null);
        bc->set_bc_name("no slip");
        bc->set_type(WhetStone::DOF_Type::POINT);
        bcs_.push_back(bc);
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

        auto bc = bc_factory.Create(spec, "kinematic", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("kinematic");
        bc->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bcs_.push_back(bc);
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

        auto bc = bc_factory.Create(spec, "traction", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("traction");
        bc->set_type(WhetStone::DOF_Type::POINT);
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

        auto bc = bc_factory.Create(spec, "normal traction", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("normal traction");
        bc->set_type(WhetStone::DOF_Type::NORMAL_COMPONENT);
        bcs_.push_back(bc);
      }
    }
  }

  // Populate matrix and preconditioner
  // -- setup phase
  double E = (*S_->Get<CV_t>(young_modulus_key_, Tags::DEFAULT).ViewComponent("cell"))[0][0];
  double nu = (*S_->Get<CV_t>(poisson_ratio_key_, Tags::DEFAULT).ViewComponent("cell"))[0][0];
  double mu =  E / (2 * (1 + nu));
  double lambda = (dim_ == 3 ) ? E * nu / (1 + nu) / (1 - 2 * nu) : E * nu / (1 + nu) / (1 - nu);

  Amanzi::WhetStone::Tensor K(dim_, 4);
  for (int i = 0; i < dim_ * dim_; ++i) K(i, i) = 2 * mu;
  for (int i = 0; i < dim_; ++i) {
    for (int j = 0; j < dim_; ++j) K(i, j) += lambda;
  }

  op_matrix_ = op_matrix_elas_->global_operator();
  op_matrix_->Init();
  op_matrix_elas_->SetTensorCoefficient(K);

  op_preconditioner_ = op_preconditioner_elas_->global_operator();
  op_preconditioner_->Init();
  op_preconditioner_elas_->SetTensorCoefficient(K);

  // -- initialize boundary conditions (memory allocation)
  auto bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::POINT));
  op_bcs_.push_back(bc);

  bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR));
  op_bcs_.push_back(bc);

  for (auto bc : op_bcs_) {
    op_matrix_elas_->AddBCs(bc, bc);
    op_preconditioner_elas_->AddBCs(bc, bc);
  }

  UpdateSourceBoundaryData_(t_ini, t_ini);
  op_matrix_elas_->EnforceBCs(*solution_);

  // -- assemble phase
  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  op_preconditioner_elas_->UpdateMatrices();
  op_preconditioner_elas_->ApplyBCs(true, true, true);

  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  op_preconditioner_elas_->global_operator()->set_inverse_parameters(pc_name,
                                                                     *preconditioner_list_);

  std::string name = ti_list_->get<std::string>("preconditioner enhancement", "none");
  if (name != "none") {
    AMANZI_ASSERT(linear_solver_list_->isSublist(name));
    Teuchos::ParameterList tmp = linear_solver_list_->sublist(name);
    op_pc_solver_ = AmanziSolvers::createIterativeMethod(tmp, op_preconditioner_);
  } else {
    op_pc_solver_ = op_preconditioner_;
  }
  op_pc_solver_->InitializeInverse();

  // we set up operators and can trigger re-initialization of stress
  Teuchos::rcp_dynamic_cast<HydrostaticStressEvaluator>(eval_hydro_stress_)->set_op(op_matrix_elas_);
  Teuchos::rcp_dynamic_cast<VolumetricStrainEvaluator>(eval_vol_strain_)->set_op(op_matrix_elas_);
  eval_->SetChanged();

  // summary of initialization
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << " TI:\"" << ti_method_name.c_str() << "\"" << std::endl
               << "matrix: "
               << op_matrix_->PrintDiagnostics() << std::endl
               << "preconditioner: "
               << op_preconditioner_->PrintDiagnostics() << std::endl;

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

  // save a copy of primary and conservative fields
  std::vector<std::string> fields({ displacement_key_, vol_strain_key_ });

  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);
  archive.CopyFieldsToPrevFields(fields, "");

  // initialization
  if (num_itrs_ == 0) {
    Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln_));
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

    // recover the original primary solution
    archive.Restore("");
    eval_->SetChanged();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Reverted displacement" << std::endl;

    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae_->CommitSolution(dt_, soln_);
  eval_->SetChanged();

  num_itrs_++;
  dt_ = dt_next_;

  return failed;
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
  if (type == Operators::OPERATOR_MATRIX)
    return op_matrix_;
  else if (type == Operators::OPERATOR_PRECONDITIONER_RAW)
    return op_preconditioner_;
  return Teuchos::null;
}

} // namespace Mechanics
} // namespace Amanzi
