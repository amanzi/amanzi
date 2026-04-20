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
#include "CommonDefs.hh"
#include "Flow_PK.hh"
#include "Energy_PK.hh"
#include "EnergyPressureEnthalpy_PK.hh"
#include "InverseFactory.hh"
#include "ModelAssumptions.hh"
#include "OperatorDefs.hh"
#include "UpwindFactory.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "PDE_DiffusionFactory.hh"
#include "StateArchive.hh"
#include "TreeOperator.hh"

#include "FlowEnergyPH_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

static constexpr double units = CommonDefs::ENTHALPY_FACTOR;

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


  compute_scaling_completed_ = false;

  left_scaling_ = my_list_->get<bool>("left scaling", false);
  left_scaling_eps_ = my_list_->get<double>("left scaling eps", 1e-4);
  right_scaling_ = my_list_->get<bool>("right scaling", false);
  P0_ = my_list_->get<double>("P0 scaling factor", 1.0);
  H0_ = my_list_->get<double>("H0 scaling factor", 1.0);

  
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
    S_->Require<CV_t, CVS_t>(state_key_, Tags::DEFAULT, state_key_, Evaluators::TS97_names)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, Evaluators::TS97_t_size);
    S_->RequireEvaluator(state_key_, Tags::DEFAULT);
  }

  // temperature
  if (!S_->HasRecord(temperature_key_)) {
    S_->Require<CV_t, CVS_t>(temperature_key_, Tags::DEFAULT, temperature_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
    S_->RequireEvaluator(temperature_key_, Tags::DEFAULT);
  }

  // densities
  if (!S_->HasRecord(mol_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    S_->RequireEvaluator(mol_density_liquid_key_, Tags::DEFAULT);
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
    S_->RequireEvaluator(iso_compressibility_key_, Tags::DEFAULT);
  }

  // viscosity
  if (!S_->HasRecord(viscosity_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(viscosity_liquid_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

    S_->RequireEvaluator(viscosity_liquid_key_, Tags::DEFAULT);
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

  std::string precon_string = my_list_->get<std::string>("preconditioner type", "full");
  //std::string precon_string = my_list_->get<std::string>("preconditioner type", "finite difference");
  
  
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "block diagonal") {
    precon_type_ = PRECON_BLOCK_DIAGONAL;
  } else if (precon_string == "full") {
    precon_type_ = PRECON_FULL;
  } else if (precon_string == "finite difference") {
    precon_type_ = PRECON_FINITEDIFF;    
  } else {
    Errors::Message message(std::string("Invalid preconditioner type ") + precon_string);
    Exceptions::amanzi_throw(message);
  }
  
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
  pde01_diff_ = opfactory.Create(oplist, mesh_);
  op01_ = pde01_diff_->global_operator();
  pde01_adv_ = opfactory_adv.Create(oplist_adv, op01_);
  op_tree_pc_->set_operator_block(0, 1, op01_);
  pde01_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op01_));

  // -- enthalpy-pressure
  pde10_diff_cond_ = opfactory.Create(oplist, mesh_);
  op10_ = pde10_diff_cond_->global_operator();
  op_tree_pc_->set_operator_block(1, 0, op10_);

  auto pks = glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  Teuchos::ParameterList upw_list = glist_->sublist("PKs").sublist(pks[0]).sublist("relative permeability");

  Operators::UpwindFactory upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, upw_list);     
  oplist.set<std::string>("nonlinear coefficient", "upwind: face");
    
  pde10_diff_flux_ = opfactory.Create(oplist, op10_);
  pde10_diff_flux_->SetTensorCoefficient(flow_pk->getK());

  pde10_adv_ = opfactory_adv.Create(oplist_adv, op10_);
  pde10_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, op10_));
 
  // solver parameters
  auto solver_list = Teuchos::sublist(glist_, "solvers");
  auto ti_list = Teuchos::sublist(my_list_, "time integrator", true);

  std::string solver_name, pc_name;
  use_cptr_prec_ = ti_list->get<bool>("use CPTR preconditioner", false);

  solver_name = ti_list->get<std::string>("linear solver", "none");
  pc_name = ti_list->get<std::string>("preconditioner", "");

  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
    pc_name, *preconditioner_list_, solver_name, *solver_list, true);
  op_tree_pc_->set_inverse_parameters(inv_list);

  if (use_cptr_prec_) {
    op_tree_amg_ = op_tree_pc_->Clone();
    op_tree_ilu_ = op_tree_pc_->Clone();

    pc_name = "ILU";
    inv_list = AmanziSolvers::mergePreconditionerSolverLists(
      pc_name, *preconditioner_list_, solver_name, *solver_list, true);
    op_tree_ilu_->set_inverse_parameters(inv_list);

    pc_name = "Hypre AMG";
    inv_list = AmanziSolvers::mergePreconditionerSolverLists(
      pc_name, *preconditioner_list_, solver_name, *solver_list, true);
    op_tree_amg_->set_inverse_parameters(inv_list);
 
    op_tree_amg_->set_block(0, 1, Teuchos::null);
    auto block = op_tree_amg_->get_operator_block(1, 1);
    /*
    for (int loop = 0; loop < 2; ++loop) {
      auto pos = block->FindMatrixOp(
        Operators::OPERATOR_SCHEMA_BASE_FACE + Operators::OPERATOR_SCHEMA_DOFS_CELL,
        Operators::OPERATOR_SCHEMA_RULE_EXACT,
        false);
      if (pos != block->end()) block->OpErase(pos);
    }
    */
  }

  // output of initialization statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (use_cptr_prec_) {
      *vo_->os() << std::endl
                 << "CPTR preconditioner:" << std::endl
                 << op_tree_amg_->PrintDiagnostics();
    }
    *vo_->os() << std::endl
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
  compute_scaling_completed_ = false;

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

  auto energy_pk = Teuchos::rcp_dynamic_cast<Energy::EnergyPressureEnthalpy_PK>(sub_pks_[1]);

  auto bc_enth = S_->GetPtrW<Operators::BCs>(bcs_enthalpy_key_, Tags::DEFAULT, "state");
  auto bc_pres = S_->GetPtrW<Operators::BCs>(bcs_flow_key_, Tags::DEFAULT, "state");

  std::string passwd("");
  Tag tag = Tags::DEFAULT;

  S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd);
  S_->GetEvaluator(viscosity_liquid_key_).Update(*S_, passwd);

  const auto& rho = S_->Get<CV_t>(mol_density_liquid_key_, tag);
  const auto& mu = S_->Get<CV_t>(viscosity_liquid_key_, tag);
  auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, tag);


  if (precon_type_ == PRECON_FULL) {
  
    // ---------------------
    // pressure-energy block
    // ---------------------
    op01_->Init();

    // -- accumulation
    S_->GetEvaluator(water_storage_key_).UpdateDerivative(*S_, passwd, enthalpy_key_, tag);
    const auto& dwsdh = S_->GetDerivative<CV_t>(water_storage_key_, tag, enthalpy_key_, tag);
    pde01_acc_->AddAccumulationTerm(dwsdh, dt, "cell");

    // -- advection div(q kt dh)
    auto coef = Teuchos::rcp(new CompositeVector(
                                                 S_->GetDerivative<CV_t>(mol_density_liquid_key_, tag, enthalpy_key_, tag)));
    coef->ReciprocalMultiply(1.0, rho, *coef, 0.0);

    S_->GetEvaluator(viscosity_liquid_key_).UpdateDerivative(*S_, passwd, enthalpy_key_, tag);
    auto coef1 = S_->GetDerivative<CV_t>(viscosity_liquid_key_, tag, enthalpy_key_, tag);
    coef->ReciprocalMultiply(-1.0, mu, coef1, 1.0);

    pde01_adv_->Setup(*flux);
    pde01_adv_->UpdateMatrices(flux.ptr(), coef.ptr());
    pde01_adv_->SetBCs(bc_enth, bc_pres);
    pde01_adv_->ApplyBCs(false, true, false);

    // ---------------------
    // energy-pressure block
    // ---------------------
    op10_->Init();

    // -- diffusion due to heat conduction
    S_->GetEvaluator(conductivity_key_).Update(*S_, passwd);
    const auto& conductivity = S_->Get<CV_t>(conductivity_key_, tag);
    coef = Teuchos::rcp(new CompositeVector(conductivity));

    S_->GetEvaluator(temperature_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
    const auto& dTdp = S_->GetDerivative<CV_t>(temperature_key_, tag, pressure_key_, tag);
    coef->Multiply(1.0, *coef, dTdp, 0.0);


    // pde10_diff_cond_->SetScalarCoefficient(coef, Teuchos::null);
    // pde10_diff_cond_->UpdateMatrices(Teuchos::null, up->Data().ptr());
    // pde10_diff_cond_->SetBCs(bc_pres, bc_enth);
    // pde10_diff_cond_->ApplyBCs(true, true, false);

  // -- diffusion due to heat transport div(K eta h / mu grad dp) (we assume relperm = 1)
    S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd);
    *coef = S_->Get<CV_t>(mol_density_liquid_key_, tag);
    coef->Multiply(1.0, *up->SubVector(1)->Data(), *coef, 0.0);
    coef->ReciprocalMultiply(1.0, mu, *coef, 0.0);


    // -- diffusion due to heat transport (we assume relperm = 1)
    // S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd);
    // *coef = S_->Get<CV_t>(mol_density_liquid_key_, tag);
    // coef->Multiply(1.0, *up->SubVector(1)->Data(), *coef, 0.0);
    // coef->ReciprocalMultiply(1.0, mu, *coef, 0.0);
  
    auto coef_wf = Teuchos::rcp(new CompositeVector(*up->SubVector(1)->Data()));
    coef_wf->Multiply(1.0, rho, *coef_wf, 0.0);
    coef_wf->ReciprocalMultiply(1.0, mu, *coef_wf, 0.0);
    *coef_wf->ViewComponent("cell") =  *coef->ViewComponent("cell");
    upwind_->Compute(*flux, bc_enth->bc_model(), *coef_wf);
  
    pde10_diff_flux_->SetScalarCoefficient(coef_wf, Teuchos::null);
    pde10_diff_flux_->UpdateMatrices(Teuchos::null, up->Data().ptr());
    pde10_diff_flux_->SetBCs(bc_pres, bc_enth);
    pde10_diff_flux_->ApplyBCs(true, true, false);
    RemoveFluxContinuityEquations_(pde10_diff_flux_);

    // -- advection div((q kt h) dp)
    S_->GetEvaluator(iso_compressibility_key_).Update(*S_, passwd);
    const auto& compressibility = S_->Get<CV_t>(iso_compressibility_key_, tag);
    coef = Teuchos::rcp(new CompositeVector(compressibility));
    coef->Multiply(1.0, *coef, *up->SubVector(1)->Data(), 0.0);

    S_->GetEvaluator(viscosity_liquid_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
    coef1 = S_->GetDerivative<CV_t>(viscosity_liquid_key_, tag, pressure_key_, tag);
    coef1.Multiply(1.0, *up->SubVector(1)->Data(), coef1, 0.0);
    coef->ReciprocalMultiply(-1.0, mu, coef1, 1.0);

    pde10_adv_->Setup(*flux);
    pde10_adv_->UpdateMatrices(flux.ptr(), coef.ptr());
    pde10_adv_->SetBCs(bc_pres, bc_enth);
    pde10_adv_->ApplyBCs(false, true, false);

    // -- accumulation
    //    modified Jacobian is used in region 4
    S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd, pressure_key_, tag);
    auto& dEdp = S_->GetDerivative<CV_t>(energy_key_, tag, pressure_key_, tag);
    pde10_acc_->AddAccumulationTerm(dEdp, dt, "cell");
  }else if (precon_type_ == PRECON_FINITEDIFF) {

    //ChangedSolution();
    //PreconditionerBlockFD_   (0, 0, t, up, dt, ddivq_dP_, pde00_adv_);
    //ChangedSolution();
    //PreconditionerAdvBlockFD_(0, 0, t, up, dt, pde00_adv_);


    // ---------------------
    // energy-pressure block
    // ---------------------
    //op10_->Init();
    
    //ChangedSolution();
    PreconditionerBlockFD_   (1, 0, t, up, dt, pde10_diff_flux_, pde10_adv_);
    //ChangedSolution();    
    PreconditionerAdvBlockFD_(1, 0, t, up, dt, pde10_adv_);


    // ---------------------
    // pressure-energy block
    // ---------------------
    //op01_->Init();
    
    //ChangedSolution();
    PreconditionerBlockFD_   (0, 1, t, up, dt, pde01_diff_, pde01_adv_);
    //ChangedSolution();
    PreconditionerAdvBlockFD_(0, 1, t, up, dt, pde01_adv_);

    //ChangedSolution();
    //PreconditionerBlockFD_   (1, 1, t, up, dt, ddivKgT_dT_, pde11_adv_);
    //ChangedSolution();
    //PreconditionerAdvBlockFD_(1, 1, t, up, dt, pde11_adv_);
    
  }

  if (use_cptr_prec_) {
    op_tree_ilu_->AssembleMatrix();
    op_tree_ilu_->InitializeInverse();
    op_tree_ilu_->ComputeInverse();  // FIXME shoud be called automatically by ApplyInverse()
  } else {
    op_tree_pc_->AssembleMatrix();
    op_tree_pc_->InitializeInverse();
    op_tree_pc_->ComputeInverse();
  }
}


/* *******************************************************************
* Selection of default or full preconditioner
******************************************************************* */
int
FlowEnergyPH_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{


  // int ierr;

  // for (int c=0;c<20; ++c) std::cout<<(*X->SubVector(1)->Data()->ViewComponent("cell"))[0][c]<<"\n";

  // ierr = op_tree_pc_->ApplyInverse(*X, *Y);

  // for (int c=0;c<20; ++c) std::cout<<(*Y->SubVector(1)->Data()->ViewComponent("cell"))[0][c]<<"\n";
  
  // if (right_scaling_) {
  //   Teuchos::RCP<TreeVector> flow_du = Y->SubVector(0);
  //   Teuchos::RCP<TreeVector> ener_du = Y->SubVector(1);

  //   flow_du->Scale(1./P0_);
  //   ener_du->Scale(1./H0_);
  // }
  
  
  // return ierr;
  // // return PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);

  int ok;
  if (!use_cptr_prec_) {
    ok = op_tree_pc_->ApplyInverse(*X, *Y);
    // return PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);
  } else {
    // Operators::Impl::TreeOperator_BlockDiagonalPreconditioner gs(*op_tree_amg_);
    Operators::Impl::TreeOperator_BlockTriangularPreconditioner gs(*op_tree_amg_);
    ok = gs.ApplyInverse(*X, *Y);

    TreeVector res(*Y), Y2(*Y);
    op_tree_ilu_->Apply(*Y, res);
    res.Update(1.0, *X, -1.0);  // r = x - J * inv(J_ell) x

    ok = op_tree_ilu_->ApplyInverse(res, Y2);
    Y->Update(1.0, Y2, 1.0);
  }
  return ok;

}


/********************************************************************
* Modifies nonlinear update du based on the maximum allowed change
* of saturation.
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
FlowEnergyPH_PK::ModifyCorrection(double dt,
                                  Teuchos::RCP<const TreeVector> f,
                                  Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> du)
{

  const auto& state_c = *S_->Get<CV_t>(state_key_).ViewComponent("cell");
  const auto& p_c = *u->SubVector(0)->Data()->ViewComponent("cell");
  const auto& h_c = *u->SubVector(1)->Data()->ViewComponent("cell");

  const auto& dp_c = *du->SubVector(0)->Data()->ViewComponent("cell");
  const auto& dh_c = *du->SubVector(1)->Data()->ViewComponent("cell");

  int nclipped = 0;
  int ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                           AmanziMesh::Parallel_kind::OWNED);

  if (right_scaling_) {
    du->SubVector(0)->Data()->Scale(P0_);
    du->SubVector(1)->Data()->Scale(H0_);
  }

  
  // increment clipping

  double max_change(0.02);
  for (int c = 0; c < ncells_owned; ++c) {
    double tmp = std::fabs(h_c[0][c]) * max_change;
    //if (c<20) std::cout<<h_c[0][c]<<" "<<dh_c[0][c]<<" "<<tmp<<"\n";

    dh_c[0][c] = std::clamp(dh_c[0][c], -tmp, tmp);


    tmp = std::fabs(p_c[0][c]) * max_change;
    //if (c<20) std::cout<<p_c[0][c]<<" "<<dp_c[0][c]<<" "<<tmp<<"\n";
    dp_c[0][c] = std::clamp(dp_c[0][c], -tmp, tmp);
    if (std::fabs(std::fabs(dp_c[0][c]) - tmp) < 1e-8 * tmp) nclipped++;
  }


  // return (nclipped) > 0 ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED :
  //                         AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;


  // phase-change checks
  double p1, h1, hf, hg, dhf, dhg, hkJ, pMPa, eps1(1e-10), eps2(1e-4);
  double units(CommonDefs::ENTHALPY_FACTOR);

  Teuchos::ParameterList plist;
  AmanziEOS::IAPWS97 eos(plist);

  for (int c = 0; c < ncells_owned; ++c) {
    p1 = p_c[0][c] - dp_c[0][c];
    h1 = h_c[0][c] - dh_c[0][c];
    //std::cout<<c<<" p1 "<<p1<<" h1 "<<h1<<"\n";
    
    pMPa = p1 * 1.0e-6;
    hkJ = h1 / units;
    auto [prop1, liquid1, vapor1] = eos.ThermodynamicsPH(pMPa, hkJ);    
    
    int rgn0 = state_c[0][c];
    int rgn1 = (int)prop1.rgn;

    //std::cout<<"rgn0 "<<rgn0<<" rgn1 "<<rgn1<<"\n";

    // leaving region 4
    if (rgn0 == 4 && rgn1 != 4) {
      pMPa = p_c[0][c] * 1.0e-6;
      hkJ = h_c[0][c] / units;
      auto [prop0, liquid0, vapor0] = eos.ThermodynamicsPH(pMPa, hkJ);

      hf = liquid0.h * units;
      hg = vapor0.h * units;
      dhf = std::max(eps1 * hg, 10.0 * eps2 * (hg - hf));
      dhg = std::max(eps1 * hg, 100.0 * eps2 * (hg - hf));

      // far outside region 4 is fine
      if (h1 <= hf - dhf || h1 >= hg + dhg) {
        nclipped++;
      } else if (rgn1 == 1 || rgn1 == 3) {
        dh_c[0][c] = h_c[0][c] - (hf + dhf);
      } else if (rgn1 == 2) {
        dh_c[0][c] = h_c[0][c] - (hg - dhg);
      } 

    // entering region 4
    } else if (rgn1 == 4 && rgn0 != 4) {
      pMPa = p_c[0][c] * 1.0e-6;
      hkJ = h_c[0][c] / units;
      auto [prop0, liquid0, vapor0] = eos.ThermodynamicsPH(pMPa, hkJ);

      hf = liquid0.h * units;
      hg = vapor0.h * units;
      dhf = std::max(eps1 * hg, eps2 * (hg - hf));
      dhg = dhf;

      // deep inside region 4 is fine
      // close to the saturation boundary kick back to the orginal region
      if (hf + dhf <= h1 && h1 <= hg - dhg) {
        nclipped++;
      } else if (rgn0 == 1 || rgn0 == 3) {
        dh_c[0][c] = h_c[0][c] - (hf - dhf);
      } else if (rgn0 == 2) {
        dh_c[0][c] = h_c[0][c] - (hg + dhg);
      } 
    }
  } 

  //  std::cout<<"nclipped "<<nclipped<<"\n";
  
  return (nclipped) > 0 ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED :
                          AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}


void FlowEnergyPH_PK::PreconditionerBlockFD_(int idi, int idj, double t,
                                           Teuchos::RCP<const TreeVector> up, double dt,
                                           Teuchos::RCP<Operators::PDE_Diffusion> pde_block,
                                           Teuchos::RCP<Operators::PDE_Advection> pde_adv){
  
    Teuchos::RCP<TreeVector> unew = Teuchos::rcp(new TreeVector(*up));
    unew->SubVector(0)->SetData(S_->GetPtrW<CompositeVector>(pressure_key_, Tags::DEFAULT, ""));
    unew->SubVector(1)->SetData(S_->GetPtrW<CompositeVector>(enthalpy_key_, Tags::DEFAULT, ""));
    
    Teuchos::RCP<TreeVector> f0 = Teuchos::rcp(new TreeVector(*unew));
    Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector(*unew));

    int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  
    double max, factor, eps(1e-8), t_new, t_old;
    up->SubVector(idj)->NormInf(&max);                                                                                                                                                                        
    factor = eps * max;
    t_old = t;
    t_new = t_old + dt;

    auto u0_tv = unew->SubVector(idj);
    auto u0_cv = u0_tv->Data();      
    auto& u0_c = *u0_cv->ViewComponent("cell");
    auto& u0_f = *u0_cv->ViewComponent("face");

    sub_pks_[0]->ChangedSolution();
    sub_pks_[1]->ChangedSolution();        
  
    FunctionalResidual(t_old, t_new, unew, unew, f0);
    
    for (int c=0; c<ncells; c++){
      //std::cout<<"cell "<<c<<"\n";
      const auto& faces = mesh_->getCellFaces(c);
      int nf = faces.size();
      WhetStone::DenseMatrix A(nf+1, nf+1);
      A.PutScalar(0.0);
      
      u0_c[0][c] += factor;
      ChangedSolution();
      
      FunctionalResidual(t_old, t_new, unew, unew, f1);
      u0_c[0][c] -= factor;
      ChangedSolution();

      auto& f0_p_c = *f0->SubVector(idi)->Data()->ViewComponent("cell");
      auto& f1_p_c = *f1->SubVector(idi)->Data()->ViewComponent("cell");
      auto& f0_p_f = *f0->SubVector(idi)->Data()->ViewComponent("face");
      auto& f1_p_f = *f1->SubVector(idi)->Data()->ViewComponent("face");

      
      A(nf,nf) = (f1_p_c[0][c] - f0_p_c[0][c]) / factor;
      
      for (int i=0; i<nf;i++){
        int f = faces[i];
        A(i,nf) = (f1_p_f[0][f] - f0_p_f[0][f]) / factor;
      }

      for (int i=0; i<nf;i++){
        sub_pks_[1]->ChangedSolution();
        int fi = faces[i];
        u0_f[0][fi] += factor;
        
        FunctionalResidual(t_old, t_new, unew, unew, f1);
        u0_f[0][fi] -= factor;

        A(nf, i) = (f1_p_c[0][c] - f0_p_c[0][c]) / factor;
        
        for (int j=0; j<nf;j++){
          int fj = faces[j];
          A(j,i) = (f1_p_f[0][fj] - f0_p_f[0][fj]) / factor;
          if (i==j){
            const auto& cells = mesh_->getFaceCells(fj); 
            if (cells.size()==2) A(j,j) *= 0.5;
          }
        }
      }


      pde_block->SetMatrix(A, c);

      // for (int i=0; i<nf; i++){
      //   int f = faces[i];
      //   const auto& cells = mesh_->getFaceCells(f);
      //   if (cells.size()==2){
      //     WhetStone::DenseMatrix A(2,2);
      //     A.PutScalar(0.0);
      //     auto u0_tv = unew->SubVector(idj);
      //     auto u0_cv = u0_tv->Data();
      
      //     auto& u0_c = *u0_cv->ViewComponent("cell");
      //     for (int i=0; i<2; i++){
      //       u0_c[0][cells[i]] += factor;
          
      //       sub_pks_[0]->ChangedSolution();
      //       sub_pks_[1]->ChangedSolution();
          
      //       FunctionalResidual(t_old, t_new, unew, unew, f1);
      //       u0_c[0][cells[i]] -= factor;

      //       sub_pks_[0]->ChangedSolution();
      //       sub_pks_[1]->ChangedSolution();

      //       auto& f0_p_c = *f0->SubVector(idi)->Data()->ViewComponent("cell");
      //       auto& f1_p_c = *f1->SubVector(idi)->Data()->ViewComponent("cell");
      //       A(1-i,i) = (f1_p_c[0][cells[1-i]] - f0_p_c[0][cells[1-i]]) / factor;          
      //     }
      //   }
      // }

    }

    //    if (idi==1&&idj==1) exit(0);
}


  void FlowEnergyPH_PK::PreconditionerAdvBlockFD_(int idi, int idj, double t,
                                           Teuchos::RCP<const TreeVector> up, double dt,
                                           Teuchos::RCP<Operators::PDE_Advection> pde_adv){
  
    Teuchos::RCP<TreeVector> unew = Teuchos::rcp(new TreeVector(*up));
    unew->SubVector(0)->SetData(S_->GetPtrW<CompositeVector>(pressure_key_, Tags::DEFAULT, ""));
    unew->SubVector(1)->SetData(S_->GetPtrW<CompositeVector>(enthalpy_key_, Tags::DEFAULT, ""));
    
    Teuchos::RCP<TreeVector> f0 = Teuchos::rcp(new TreeVector(*unew));
    Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector(*unew));

    int nfaces = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  
    double max, factor, eps(1e-8), t_new, t_old;
    up->SubVector(idj)->NormInf(&max);                                                                                                                                                                        
    factor = eps * max;
    t_old = t;
    t_new = t_old + dt;

    sub_pks_[0]->ChangedSolution();
    sub_pks_[1]->ChangedSolution();
  
    FunctionalResidual(t_old, t_new, unew, unew, f0);

    for (int f=0; f<nfaces; f++){
      const auto& cells = mesh_->getFaceCells(f);
      if (cells.size()==2){
        WhetStone::DenseMatrix A(2,2);
        A.PutScalar(0.0);
        auto u0_tv = unew->SubVector(idj);
        auto u0_cv = u0_tv->Data();
      
        auto& u0_c = *u0_cv->ViewComponent("cell");
        for (int i=0; i<2; i++){
          u0_c[0][cells[i]] += factor;

          ChangedSolution();          
          FunctionalResidual(t_old, t_new, unew, unew, f1);

          auto& f0_p_c = *f0->SubVector(idi)->Data()->ViewComponent("cell");
          auto& f1_p_c = *f1->SubVector(idi)->Data()->ViewComponent("cell");          
          
          A(1-i,i) = (f1_p_c[0][cells[1-i]] - f0_p_c[0][cells[1-i]]) / factor;
          
          u0_c[0][cells[i]] -= factor;          
          ChangedSolution();
          
        }
        
        pde_adv->SetMatrix(A, f);
        
      }      
    }
  
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
