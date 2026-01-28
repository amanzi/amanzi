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

#include "BCs.hh"
#include "Flow_PK.hh"
#include "Energy_PK.hh"
#include "ModelAssumptions.hh"
#include "OperatorDefs.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "PDE_DiffusionFactory.hh"
#include "StateArchive.hh"

#include "FlowEnergyPH_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
FlowEnergyPH_PK::FlowEnergyPH_PK(Teuchos::ParameterList& pk_tree,
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
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  domain_ = my_list_->template get<std::string>("domain name", "domain");

  vo_ = Teuchos::rcp(new VerboseObject("FlowEnergy-" + domain_, *my_list_));
}


/* *******************************************************************
* Physics-based setup of PK.
******************************************************************* */
void
FlowEnergyPH_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);

  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(my_list_, "physical models and assumptions");
  ModelAssumptions assumptions;
  assumptions.Init(*physical_models, *mesh_);

  // keys
  pressure_key_ = Keys::getKey(domain_, "pressure");
  enthalpy_key_ = Keys::getKey(domain_, "enthalpy");
  temperature_key_ = Keys::getKey(domain_, "temperature");
  energy_key_ = Keys::getKey(domain_, "energy");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");

  state_key_ = Keys::getKey(domain_, "thermodynamic_state");
  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid");
  iso_compressibility_key_ = Keys::getKey(domain_, "isothermal_compressibility");
  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");

  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");
  water_storage_key_ = Keys::getKey(domain_, "water_storage");
  bcs_flow_key_ = Keys::getKey(domain_, "bcs_flow");
  bcs_enthalpy_key_ = Keys::getKey(domain_, "bcs_enthalpy");

  // rock
  if (!S_->HasRecord(particle_density_key_)) {
    S_->Require<CV_t, CVS_t>(particle_density_key_, Tags::DEFAULT, particle_density_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(particle_density_key_, Tags::DEFAULT);
  }

  // thermodynamics (need a few evaluators here to enforce IAPWS97 formulation)
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

  // temperature
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

  // densities
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

    S_->RequireDerivative<CV_t, CVS_t>(mol_density_liquid_key_,
                                       Tags::DEFAULT,
                                       enthalpy_key_,
                                       Tags::DEFAULT,
                                       mol_density_liquid_key_).SetGhosted();
  }

  // isothermal compressibility
  if (!S_->HasRecord(iso_compressibility_key_)) {
    S_->Require<CV_t, CVS_t>(iso_compressibility_key_, Tags::DEFAULT, iso_compressibility_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    Teuchos::ParameterList elist(iso_compressibility_key_);
    elist.set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new Evaluators::IsothermalCompressibilityEvaluator(elist));
    S_->SetEvaluator(iso_compressibility_key_, Tags::DEFAULT, eval);
  }

  // viscosity
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
        viscosity_liquid_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetGhosted();

    S_->RequireDerivative<CV_t, CVS_t>(
        viscosity_liquid_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetGhosted();
  }

  // inform other PKs about strong coupling
  // -- flow
  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  Teuchos::ParameterList& flow =
    glist_->sublist("PKs").sublist(pks[0]).sublist("physical models and assumptions");
  flow.set<bool>("thermoelasticity", assumptions.thermoelasticity)
    .set<bool>("biot scheme: undrained split", assumptions.undrained_split)
    .set<bool>("biot scheme: fixed stress split", assumptions.fixed_stress_split);

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup();

  // extend state structure
  S_->RequireDerivative<CV_t, CVS_t>(
      water_storage_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT, water_storage_key_)
    .SetGhosted();
}


/* *******************************************************************
* Initialization of copies requires fileds to exists
******************************************************************* */
void
FlowEnergyPH_PK::Initialize()
{
  auto flow_pk = Teuchos::rcp_dynamic_cast<Flow::Flow_PK>(sub_pks_[0]);

  // complete initialization of variables
  Amanzi::DeriveFaceValuesFromCellValues(S_->GetW<CV_t>(viscosity_liquid_key_, Tags::DEFAULT, viscosity_liquid_key_));
  Amanzi::DeriveFaceValuesFromCellValues(S_->GetW<CV_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_));

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

  // cross coupling
  Operators::PDE_DiffusionFactory opfactory;
  Teuchos::ParameterList oplist = Teuchos::rcp_dynamic_cast<Energy::Energy_PK>(sub_pks_[1])
    ->getPList()->sublist("operators").sublist("diffusion operator").sublist("preconditioner");

  Operators::PDE_AdvectionUpwindFactory opfactory_adv;
  Teuchos::ParameterList oplist_adv = Teuchos::rcp_dynamic_cast<Energy::Energy_PK>(sub_pks_[1])
    ->getPList()->sublist("operators").sublist("advection operator");

  // -- pressure-enthalpy 
  pde01_adv_ = opfactory_adv.Create(oplist_adv, mesh_);
  op01_ = pde01_adv_->global_operator();
  op_tree_pc_->set_operator_block(0, 1, op01_);
  pde01_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op01_));

  // -- enthalpy-pressure
  pde10_diff_cond_ = opfactory.Create(oplist, mesh_);
  op10_ = pde10_diff_cond_->global_operator();
  op_tree_pc_->set_operator_block(1, 0, op10_);

  pde10_diff_flux_ = opfactory.Create(oplist, op10_);
  pde10_diff_flux_->SetTensorCoefficient(flow_pk->getK());

  pde10_adv_ = opfactory_adv.Create(oplist_adv, op10_);
  pde10_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op10_));

  std::string pc_name = my_list_->sublist("time integrator").get<std::string>("preconditioner");
  op_tree_pc_->set_inverse_parameters(pc_name, *preconditioner_list_);

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
* Performs one timestep.
******************************************************************* */
bool
FlowEnergyPH_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // save a copy of conservative fields
  std::vector<std::string> fields( { energy_key_, pressure_key_, enthalpy_key_ });

  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

  // try timestep
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);
  if (fail) archive.Restore("");

  return fail;
}


/* *******************************************************************
* Performs one timestep.
******************************************************************* */
void
FlowEnergyPH_PK::FunctionalResidual(double t_old,
                                    double t_new,
                                    Teuchos::RCP<const TreeVector> u_old,
                                    Teuchos::RCP<TreeVector> u_new,
                                    Teuchos::RCP<TreeVector> f)
{
  // update molar flux
  auto u_new0 = u_new->SubVector(0);
  auto flux = S_->GetPtrW<CV_t>(mol_flowrate_key_, Tags::DEFAULT, "");
  auto op0 = sub_pks_[0]->my_pde(Operators::PDEType::PDE_DIFFUSION);
  op0->UpdateFlux(u_new0->Data().ptr(), flux.ptr());

  // flow
  auto u_old0 = u_old->SubVector(0);
  auto f0 = f->SubVector(0);
  sub_pks_[0]->FunctionalResidual(t_old, t_new, u_old0, u_new0, f0);

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
FlowEnergyPH_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, up, dt);

  auto bc_enth = S_->GetPtrW<Operators::BCs>(bcs_enthalpy_key_, Tags::DEFAULT, "state");
  auto bc_pres = S_->GetPtrW<Operators::BCs>(bcs_flow_key_, Tags::DEFAULT, "state");

  std::string passwd("");
  Tag tag = Tags::DEFAULT;

  S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd);
  S_->GetEvaluator(viscosity_liquid_key_).Update(*S_, passwd);

  const auto& rho = S_->Get<CV_t>(mol_density_liquid_key_, Tags::DEFAULT);
  const auto& mu = S_->Get<CV_t>(viscosity_liquid_key_, Tags::DEFAULT);
  auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, Tags::DEFAULT);

  // pressure-energy block
  op01_->Init();

  // -- accumulation
  S_->GetEvaluator(water_storage_key_).UpdateDerivative(*S_, passwd, enthalpy_key_, tag);
  const auto& dwsdh = S_->GetDerivative<CV_t>(water_storage_key_, tag, enthalpy_key_, tag);
  pde01_acc_->AddAccumulationTerm(dwsdh, dt, "cell");

  // -- advection div(q kt dh)
  auto coef = Teuchos::rcp(new CompositeVector(
    S_->GetDerivative<CV_t>(mol_density_liquid_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT)));
  coef->ReciprocalMultiply(1.0, rho, *coef, 0.0);

  S_->GetEvaluator(viscosity_liquid_key_).UpdateDerivative(*S_, passwd, enthalpy_key_, Tags::DEFAULT);
  auto coef1 = S_->GetDerivative<CV_t>(viscosity_liquid_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT);
  coef->ReciprocalMultiply(-1.0, mu, coef1, 1.0);

  pde01_adv_->Setup(*flux);
  pde01_adv_->UpdateMatrices(flux.ptr(), coef.ptr());
  pde01_adv_->SetBCs(bc_enth, bc_pres);
  pde01_adv_->ApplyBCs(false, true, false);

  // energy-pressure block
  op10_->Init();

  // -- diffusion due to heat conduction
  S_->GetEvaluator(conductivity_key_).Update(*S_, passwd);
  const auto& conductivity = S_->Get<CV_t>(conductivity_key_, tag);
  coef = Teuchos::rcp(new CompositeVector(conductivity));

  S_->GetEvaluator(temperature_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
  const auto& dTdp = S_->GetDerivative<CV_t>(temperature_key_, tag, pressure_key_, tag);
  coef->Multiply(1.0, *coef, dTdp, 0.0);

  pde10_diff_cond_->SetScalarCoefficient(coef, Teuchos::null);
  pde10_diff_cond_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  pde10_diff_cond_->SetBCs(bc_pres, bc_enth);
  pde10_diff_cond_->ApplyBCs(true, true, false);

  // -- diffusion due to heat transport (we assume relperm = 1)
  S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd);
  *coef = S_->Get<CV_t>(mol_density_liquid_key_, tag);
  coef->Multiply(1.0, *up->SubVector(1)->Data(), *coef, 0.0);
  coef->ReciprocalMultiply(1.0, mu, *coef, 0.0);

  pde10_diff_flux_->SetScalarCoefficient(coef, Teuchos::null);
  pde10_diff_flux_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  pde10_diff_flux_->SetBCs(bc_pres, bc_enth);
  pde10_diff_flux_->ApplyBCs(true, true, false);
  RemoveFluxContinuityEquations_(pde10_diff_flux_);

  // -- advection div((q kt h) dp)
  S_->GetEvaluator(iso_compressibility_key_).Update(*S_, passwd);
  const auto& compressibility = S_->Get<CV_t>(iso_compressibility_key_, Tags::DEFAULT);
  coef = Teuchos::rcp(new CompositeVector(compressibility));
  coef->Multiply(1.0, *coef, *up->SubVector(1)->Data(), 0.0);

  S_->GetEvaluator(viscosity_liquid_key_).UpdateDerivative(*S_, passwd, pressure_key_, Tags::DEFAULT);
  coef1 = S_->GetDerivative<CV_t>(viscosity_liquid_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT);
  coef1.Multiply(1.0, *up->SubVector(1)->Data(), coef1, 0.0);
  coef->ReciprocalMultiply(-1.0, mu, coef1, 1.0);

  pde10_adv_->Setup(*flux);
  pde10_adv_->UpdateMatrices(flux.ptr(), coef.ptr());
  pde10_adv_->SetBCs(bc_pres, bc_enth);
  pde10_adv_->ApplyBCs(false, true, false);

  // -- accumulation
  S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
  const auto& dEdp = S_->GetDerivative<CV_t>(energy_key_, tag, pressure_key_, tag);
  pde10_acc_->AddAccumulationTerm(dEdp, dt, "cell");

  if (!symbolic_assembly_complete_) {
    op_tree_pc_->SymbolicAssembleMatrix();
    symbolic_assembly_complete_ = true;
  }
  op_tree_pc_->AssembleMatrix();
  // std::cout << *op_tree_pc_->A() << std::endl; exit(0);
  op_tree_pc_->InitializeInverse();
  op_tree_pc_->ComputeInverse();
}


/* *******************************************************************
* Selection of default or full preconditioner
******************************************************************* */
int
FlowEnergyPH_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  return op_tree_pc_->ApplyInverse(*X, *Y);
  // return PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);
}


/* *******************************************************************
* L-scheme stability constant is updated by PKs.
******************************************************************* */
std::vector<Key>
FlowEnergyPH_PK::SetupLSchemeKey(Teuchos::ParameterList& plist)
{
  auto tmp = sub_pks_[0]->SetupLSchemeKey(plist);
  L_scheme_keys_.insert(L_scheme_keys_.end(), tmp.begin(), tmp.end());
  return L_scheme_keys_;
}


/* *******************************************************************
* Selection of default or full preconditioner
******************************************************************* */
void FlowEnergyPH_PK::RemoveFluxContinuityEquations_(Teuchos::RCP<Operators::PDE_Diffusion>& pde)
{
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    WhetStone::DenseMatrix& Acell = pde->local_op()->matrices[c];

    int nfaces = mesh_->getCellNumFaces(c);
    for (int m = 0; m < nfaces; ++m) {
      for (int n = 0; n < nfaces + 1; ++n) {
        Acell(m, n) = 0.0;
      }
    }
  } 
}

} // namespace Amanzi
