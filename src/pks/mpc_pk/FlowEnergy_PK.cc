/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling Flow PK with Energy PK.
*/

#include "Flow_PK.hh"
#include "Energy_PK.hh"
#include "OperatorDefs.hh"
#include "StateArchive.hh"

#include "FlowEnergy_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
FlowEnergy_PK::FlowEnergy_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln)
  : Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln),
    glist_(glist)
{
  // we will use a few parameter lists
  auto pk_list = Teuchos::sublist(glist, "PKs", true);
  my_list_ = Teuchos::sublist(pk_list, name_, true);
  domain_ = my_list_->template get<std::string>("domain name", "domain");

  vo_ = Teuchos::rcp(new VerboseObject("FlowEnergy-" + domain_, *my_list_));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
FlowEnergy_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);

  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(my_list_, "physical models and assumptions");
  bool vapor_diff = physical_models->get<bool>("vapor diffusion");
  bool thermoelasticity = physical_models->get<bool>("thermoelasticity", false);

  if (physical_models->isParameter("eos lookup table")) {
    eos_table_ = physical_models->get<std::string>("eos lookup table");
  }

  // keys
  temperature_key_ = Keys::getKey(domain_, "temperature");
  energy_key_ = Keys::getKey(domain_, "energy");
  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");

  pressure_key_ = Keys::getKey(domain_, "pressure");
  sat_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  ws_key_ = Keys::getKey(domain_, "water_storage");

  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");

  // Fields for solids
  // -- rock
  if (!S_->HasRecord(particle_density_key_)) {
    S_->Require<CV_t, CVS_t>(particle_density_key_, Tags::DEFAULT, particle_density_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(particle_density_key_, Tags::DEFAULT);
  }

  // Fields for liquid
  // -- internal energy
  if (!S_->HasRecord(ie_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(ie_liquid_key_, Tags::DEFAULT, ie_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
    S_->RequireEvaluator(ie_liquid_key_, Tags::DEFAULT);
  }

  // -- molar and mass density
  if (!S_->HasRecord(mol_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
    S_->RequireEvaluator(mol_density_liquid_key_, Tags::DEFAULT);
  }

  S_->RequireEvaluator(mass_density_liquid_key_, Tags::DEFAULT);

  // inform other PKs about strong coupling
  // -- flow
  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  std::string model = (vapor_diff) ? "two-phase" : "one-phase";
  Teuchos::ParameterList& flow =
    glist_->sublist("PKs").sublist(pks[0]).sublist("physical models and assumptions");
  flow.set<bool>("vapor diffusion", vapor_diff)
    .set<bool>("thermoelasticity", thermoelasticity)
    .set<bool>("biot scheme: undrained split",
               physical_models->get<bool>("biot scheme: undrained split", false))
    .set<bool>("biot scheme: fixed stress split",
               physical_models->get<bool>("biot scheme: fixed stress split", false));

  // -- energy
  Teuchos::ParameterList& energy =
    glist_->sublist("PKs").sublist(pks[1]).sublist("physical models and assumptions");
  energy.set<bool>("vapor diffusion", vapor_diff);

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup();

  // extend state structure
  S_->RequireDerivative<CV_t, CVS_t>(
      energy_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, energy_key_)
    .SetGhosted();

  S_->RequireDerivative<CV_t, CVS_t>(
      ws_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ws_key_)
    .SetGhosted();
}


/* *******************************************************************
* Initialization of copies requires fileds to exists
******************************************************************* */
void
FlowEnergy_PK::Initialize()
{
  include_pt_coupling_ =
    my_list_->sublist("time integrator").get<bool>("include coupling terms", false);

  Amanzi::PK_MPCStrong<PK_BDF>::Initialize();

  // MPC_PKs that build on top of this may need a tree operator. Since
  // we cannot use solution_, we create a TVS from scratch
  auto op0 = sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX);
  auto op1 = sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX);

  auto tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(CreateTVSwithOneLeaf(op0->DomainMap()));
  tvs->PushBack(CreateTVSwithOneLeaf(op1->DomainMap()));

  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_matrix_->set_operator_block(0, 0, op0);
  op_tree_matrix_->set_operator_block(1, 1, op1);

  op_tree_pc_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_pc_->set_operator_block(
    0, 0, sub_pks_[0]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW));
  op_tree_pc_->set_operator_block(
    1, 1, sub_pks_[1]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW));

  op_tree_rhs_ = Teuchos::rcp(new TreeVector());
  op_tree_rhs_->PushBack(CreateTVwithOneLeaf(op0->rhs()));
  op_tree_rhs_->PushBack(CreateTVwithOneLeaf(op1->rhs()));

  // full preconditioner
  if (include_pt_coupling_) {
    op10_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
    op_tree_pc_->set_operator_block(1, 0, op10_acc_->global_operator());

    op01_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
    op_tree_pc_->set_operator_block(0, 1, op01_acc_->global_operator());

    std::string pc_name =
      Teuchos::sublist(my_list_, "time integrator", true)->get<std::string>("preconditioner");
    auto& pc_list = Teuchos::sublist(glist_, "preconditioners", true)->sublist(pc_name);
    op_tree_pc_->set_inverse_parameters(pc_list);

    op_tree_pc_->SymbolicAssembleMatrix();
  }

  // output of initialization statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl
               << "matrix:" << std::endl
               << op_tree_matrix_->PrintDiagnostics() << std::endl
               << "preconditioner:" << std::endl
               << op_tree_pc_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl
               << std::endl;
  }

  AMANZI_ASSERT(sub_pks_.size() == 2);
}


/* *******************************************************************
* Performs one time step.
******************************************************************* */
bool
FlowEnergy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // save a copy of conservative fields
  std::vector<std::string> fields(
    { pressure_key_, temperature_key_, sat_liquid_key_, energy_key_ });
  if (S_->HasRecord(ws_key_)) fields.push_back(ws_key_);

  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

  // try time step
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);
  if (fail) archive.Restore("");

  return fail;
}


/* *******************************************************************
* Performs one time step.
******************************************************************* */
void
FlowEnergy_PK::FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f)
{
  // flow
  auto u_old0 = u_old->SubVector(0);
  auto u_new0 = u_new->SubVector(0);
  auto f0 = f->SubVector(0);
  sub_pks_[0]->FunctionalResidual(t_old, t_new, u_old0, u_new0, f0);

  // update molar flux
  Key key = Keys::getKey(domain_, "molar_flow_rate");
  auto mol_flowrate = S_->GetPtrW<CV_t>(key, Tags::DEFAULT, "");
  auto op0 = sub_pks_[0]->my_pde(Operators::PDEType::PDE_DIFFUSION);
  op0->UpdateFlux(u_new0->Data().ptr(), mol_flowrate.ptr());

  if (Keys::getVarName(sub_pks_[0]->name()) == "darcy") {
    double molar_mass = S_->Get<double>("const_fluid_molar_mass");
    mol_flowrate->Scale(1.0 / molar_mass);
  }

  // energy
  auto u_old1 = u_old->SubVector(1);
  auto u_new1 = u_new->SubVector(1);
  auto f1 = f->SubVector(1);
  sub_pks_[1]->FunctionalResidual(t_old, t_new, u_old1, u_new1, f1);
}


/* *******************************************************************
* Preconditioner update
******************************************************************* */
void
FlowEnergy_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, up, dt);

  if (include_pt_coupling_) {
    std::string passwd("");
    S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd, pressure_key_, Tags::DEFAULT);
    const auto& dEdP =
      S_->GetDerivative<CV_t>(energy_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT);

    S_->GetEvaluator(ws_key_).UpdateDerivative(*S_, passwd, temperature_key_, Tags::DEFAULT);
    const auto& dws_dT =
      S_->GetDerivative<CV_t>(ws_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);

    if (dt > 0.0) {
      op10_acc_->AddAccumulationDelta(*up->SubVector(0)->Data(), dEdP, dEdP, dt, "cell");
      op01_acc_->AddAccumulationDelta(*up->SubVector(1)->Data(), dws_dT, dws_dT, dt, "cell");
    }

    op_tree_pc_->AssembleMatrix();
    op_tree_pc_->InitializeInverse();
    op_tree_pc_->ComputeInverse();
  }
}


/* *******************************************************************
* Selection of default or full preconditioner
******************************************************************* */
int
FlowEnergy_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  if (include_pt_coupling_) {
    Y->PutScalar(0.0);
    return op_tree_pc_->ApplyInverse(*X, *Y);
  } else {
    return PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);
  }
}


} // namespace Amanzi
