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
#include "ModelAssumptions.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "StateArchive.hh"

#include "FlowEnergyPT_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* *******************************************************************
* Constructor
******************************************************************* */
FlowEnergyPT_PK::FlowEnergyPT_PK(Teuchos::ParameterList& pk_tree,
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
FlowEnergyPT_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);

  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(my_list_, "physical models and assumptions");
  ModelAssumptions assumptions;
  assumptions.Init(*physical_models, *mesh_);

  // keys
  temperature_key_ = Keys::getKey(domain_, "temperature");
  energy_key_ = Keys::getKey(domain_, "energy");
  enthalpy_key_ = Keys::getKey(domain_, "enthalpy");
  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");

  pressure_key_ = Keys::getKey(domain_, "pressure");
  sat_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  ws_key_ = Keys::getKey(domain_, "water_storage");

  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");
  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid");

  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");
  permeability_key_ = Keys::getKey(domain_, "permeability");
  
  aperture_key_ = Keys::getKey(domain_, "aperture");
  conductivity_eff_key_ = Keys::getKey(domain_, "thermal_conductivity_effective");
  conductivity_gen_key_ =
    (!assumptions.flow_on_manifold) ? conductivity_key_ : conductivity_eff_key_;

  
  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");
  bcs_flow_key_ = Keys::getKey(domain_, "bcs_flow");
  bcs_temperature_key_ = Keys::getKey(domain_, "bcs_temperature");
  bcs_enthalpy_key_ = Keys::getKey(domain_, "bcs_enthalpy");

  alpha_key_ = Keys::getKey(domain_, "alpha_coef");
  beta_key_ = Keys::getKey(domain_, "beta_coef");
  
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
  std::string model = (assumptions.vapor_diffusion) ? "two-phase" : "one-phase";
  Teuchos::ParameterList& flow =
    glist_->sublist("PKs").sublist(pks[0]).sublist("physical models and assumptions");
  flow.set<bool>("vapor diffusion", assumptions.vapor_diffusion)
    .set<bool>("thermoelasticity", assumptions.thermoelasticity)
    .set<bool>("biot scheme: undrained split", assumptions.undrained_split)
    .set<bool>("biot scheme: fixed stress split", assumptions.fixed_stress_split);

  // -- energy
  Teuchos::ParameterList& energy =
    glist_->sublist("PKs").sublist(pks[1]).sublist("physical models and assumptions");
  energy.set<bool>("vapor diffusion", assumptions.vapor_diffusion);

  // process other PKs
  PK_MPCStrong<PK_BDF>::Setup();

  // extend state structure
  S_->RequireDerivative<CV_t, CVS_t>(
      energy_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, energy_key_)
    .SetGhosted();

  S_->RequireDerivative<CV_t, CVS_t>(
      ws_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ws_key_)
    .SetGhosted();

  S_->RequireDerivative<CV_t, CVS_t>(
      enthalpy_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, enthalpy_key_)
    .SetGhosted();

  S_->RequireDerivative<CV_t, CVS_t>(
      beta_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, beta_key_)
    .SetGhosted();
  
}


/* *******************************************************************
* Initialization of copies requires fileds to exists
******************************************************************* */
void
FlowEnergyPT_PK::Initialize()
{
  auto flow_pk = Teuchos::rcp_dynamic_cast<Flow::Flow_PK>(sub_pks_[0]);

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
    Operators::PDE_DiffusionFactory opfactory;
    Teuchos::ParameterList oplist = Teuchos::rcp_dynamic_cast<Energy::Energy_PK>(sub_pks_[1])
      ->getPList()->sublist("operators").sublist("diffusion operator").sublist("preconditioner");

    Operators::PDE_AdvectionUpwindFactory opfactory_adv;
    Teuchos::ParameterList oplist_adv = Teuchos::rcp_dynamic_cast<Energy::Energy_PK>(sub_pks_[1])
      ->getPList()->sublist("operators").sublist("advection operator");

    // -- pressure-temperature
    pde01_adv_ = opfactory_adv.Create(oplist_adv, mesh_);
    op01_ = pde01_adv_->global_operator();
    op_tree_pc_->set_operator_block(0, 1, op01_);
    pde01_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op01_));

    // -- temperature-pressure
    pde10_diff_flux_ = opfactory.Create(oplist, mesh_);
    pde10_diff_flux_->SetTensorCoefficient(flow_pk->getK());
    op10_ = pde10_diff_flux_->global_operator();
    op_tree_pc_->set_operator_block(1, 0, op10_);

    pde10_adv_ = opfactory_adv.Create(oplist_adv, op10_);
    pde10_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op10_));

  }

  std::string pc_name =
      Teuchos::sublist(my_list_, "time integrator", true)->get<std::string>("preconditioner");
    auto& pc_list = Teuchos::sublist(glist_, "preconditioners", true)->sublist(pc_name);
    op_tree_pc_->set_inverse_parameters(pc_list);

    op_tree_pc_->SymbolicAssembleMatrix();
  

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
FlowEnergyPT_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // save a copy of conservative fields
  std::vector<std::string> fields(
    { pressure_key_, temperature_key_, sat_liquid_key_, energy_key_ });
  if (S_->HasRecord(ws_key_) ) fields.push_back(ws_key_);

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
FlowEnergyPT_PK::FunctionalResidual(double t_old,
                                    double t_new,
                                    Teuchos::RCP<const TreeVector> u_old,
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
FlowEnergyPT_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, up, dt);


  std::string passwd("");
  Tag tag = Tags::DEFAULT;

  if (include_pt_coupling_) {
    const auto& rho = S_->Get<CV_t>(mol_density_liquid_key_, tag);
    const auto& mu = S_->Get<CV_t>(viscosity_liquid_key_, tag);
    const auto& enth = S_->Get<CV_t>(enthalpy_key_, tag);
    auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, tag);

    auto bc_pres = S_->GetPtrW<Operators::BCs>(bcs_flow_key_, tag, "state");
    auto bc_temp = S_->GetPtrW<Operators::BCs>(bcs_temperature_key_, tag, "state");

  
    S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd);
    S_->GetEvaluator(enthalpy_key_).Update(*S_, passwd);

    // pressure-energy block
    op01_->Init();

    // -- accumulation
    S_->GetEvaluator(ws_key_).UpdateDerivative(*S_, passwd, temperature_key_, tag);
    const auto& dwsdT = S_->GetDerivative<CV_t>(ws_key_, tag, temperature_key_, tag);
    pde01_acc_->AddAccumulationTerm(dwsdT, dt, "cell");

    // -- advection div(q kt dT)
    S_->GetEvaluator(mol_density_liquid_key_).UpdateDerivative(*S_, passwd, temperature_key_, tag);
    auto coef = Teuchos::rcp(new CompositeVector(
                                                 S_->GetDerivative<CV_t>(mol_density_liquid_key_, tag, temperature_key_, tag)));
    coef->ReciprocalMultiply(1.0, rho, *coef, 0.0);

    S_->GetEvaluator(viscosity_liquid_key_).UpdateDerivative(*S_, passwd, temperature_key_, tag);
    auto coef1 = S_->GetDerivative<CV_t>(viscosity_liquid_key_, tag, temperature_key_, tag);
    coef->ReciprocalMultiply(-1.0, mu, coef1, 1.0);

    pde01_adv_->Setup(*flux);
    pde01_adv_->UpdateMatrices(flux.ptr(), coef.ptr());
    pde01_adv_->SetBCs(bc_temp, bc_pres);
    pde01_adv_->ApplyBCs(false, true, false);

    // energy-pressure block
    op10_->Init();

    // -- diffusion due to heat transport (with assumption relperm = 1), div([K beta] grad dp)    
    // coef->Multiply(1.0, enth, *coef, 0.0);
    // coef->ReciprocalMultiply(1.0, mu, *coef, 0.0);
    S_->GetEvaluator(beta_key_).Update(*S_, passwd);
    Teuchos::RCP<CompositeVector> coef2 = Teuchos::rcp(new CompositeVector(S_->Get<CV_t>(beta_key_, tag)));    

    pde10_diff_flux_->SetScalarCoefficient(coef2, Teuchos::null);
    pde10_diff_flux_->UpdateMatrices(Teuchos::null, up->Data().ptr());
    pde10_diff_flux_->SetBCs(bc_pres, bc_temp);
    pde10_diff_flux_->ApplyBCs(true, true, false);
    

    // // -- advection due to heat transport, div([q dH/dp] dp) FIXME we need d(qH)/dp
    S_->GetEvaluator(beta_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
    Teuchos::RCP<CompositeVector> coef3 = Teuchos::rcp(new CompositeVector(S_->GetDerivative<CV_t>(beta_key_, tag, pressure_key_, tag)));

    coef3->Multiply(1.0, mu, *coef3, 0.0);
    coef3->ReciprocalMultiply(1.0, rho, *coef3, 0.0);            
    
    pde10_adv_->Setup(*flux);
    pde10_adv_->UpdateMatrices(flux.ptr(), coef3.ptr());
    pde10_adv_->SetBCs(bc_pres, bc_temp);
    pde10_adv_->ApplyBCs(false, true, false);

    // -- accumulation
    S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
    auto& dEdp = S_->GetDerivative<CV_t>(energy_key_, tag, pressure_key_, tag);
    pde10_acc_->AddAccumulationTerm(dEdp, dt, "cell");
  }
  
  op_tree_pc_->AssembleMatrix();
  op_tree_pc_->InitializeInverse();
  op_tree_pc_->ComputeInverse();
    //}
}


/* *******************************************************************
* Selection of default or full preconditioner
******************************************************************* */
int
FlowEnergyPT_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  //  if (include_pt_coupling_) {
    Y->PutScalar(0.0);
    return op_tree_pc_->ApplyInverse(*X, *Y);
  // } else {
  //   return PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);
  // }
}

// -----------------------------------------------------------------------------
// L-scheme stability constant is updated by PKs.
// -----------------------------------------------------------------------------
std::vector<Key>
FlowEnergyPT_PK::SetupLSchemeKey(Teuchos::ParameterList& plist)
{
  auto tmp = sub_pks_[0]->SetupLSchemeKey(plist);
  L_scheme_keys_.insert(L_scheme_keys_.end(), tmp.begin(), tmp.end());
  return L_scheme_keys_;
}

} // namespace Amanzi
