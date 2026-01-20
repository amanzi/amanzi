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
#include "TreeOperator.hh"
#include "StateArchive.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_Advection.hh"
#include "PDE_Accumulation.hh"
#include "Operator.hh"
#include "FlowEnergy_PK.hh"
#include "PK_MPCStrong.hh"
#include "IO.hh"

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
  //  auto pk_order = my_list_->get<Teuchos::Array<std::string>>("PKs order");
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *my_list_));
  
  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(my_list_, "physical models and assumptions");
  ModelAssumptions assumptions;
  assumptions.Init(*physical_models, *mesh_);

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

  enth_key_ = Keys::getKey(domain_, "enthalpy");
  hkr_key_  = Keys::getKey(domain_, "enthalpy_times_relative_permeability");
  kr_key_   = Keys::getKey(domain_, "relative_permeability");
  
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


  include_pt_coupling_ =
    my_list_->sublist("time integrator").get<bool>("include coupling terms", false);

  std::string precon_string = my_list_->get<std::string>("preconditioner type", "none");
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

  // full preconditioner
  if (precon_type_ == PRECON_FULL) {

    std::vector<AmanziMesh::Entity_kind> locations2(2);
    std::vector<std::string> names2(2);
    std::vector<int> num_dofs2(2, 1);
    locations2[0] = AmanziMesh::Entity_kind::CELL;
    names2[0] = "cell";

    // need a field, evaluator, and upwinding for h * kr * rho/mu
    // -- first the evaluator
    auto& hkr_eval_list = S_->GetEvaluatorList(hkr_key_);
    hkr_eval_list.set("evaluator name", hkr_key_);
    Teuchos::Array<std::string> deps(2);
    deps[0] = enth_key_;
    deps[1] = kr_key_;
    hkr_eval_list.set("multiplicative dependency keys", deps);
    hkr_eval_list.set("evaluator type", "multiplicative reciprocal");


    
    // -- now the field
    names2[1] = "boundary_face";
    locations2[1] = AmanziMesh::Entity_kind::BOUNDARY_FACE;
    S_->Require<CompositeVector, CompositeVectorSpace>(hkr_key_, Tags::DEFAULT)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
    S_->RequireEvaluator(hkr_key_, Tags::DEFAULT);

    // S_->Require<CompositeVector, CompositeVectorSpace>(uw_hkr_key_, Tags::DEFAULT, name_)
    //   .SetMesh(mesh_)
    //   ->SetGhosted()
    //   ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
    // S_->GetRecordW(uw_hkr_key_, Tags::DEFAULT, name_).set_io_vis(false);

    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
                       hkr_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);
    S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
                       hkr_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT);     

  }
  
}


/* *******************************************************************
* Initialization of copies requires fileds to exists
******************************************************************* */
void
FlowEnergy_PK::Initialize()
{

  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  
  Amanzi::PK_MPCStrong<PK_BDF>::Initialize();

  // MPC_PKs that build on top of this may need a tree operator. Since
  // we cannot use solution_, we create a TVS from scratch
  auto op0 = sub_pks_[0]->my_operator(Operators::OPERATOR_MATRIX);
  auto op1 = sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX);

  // Get the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Operator> pc_block00 = sub_pks_[0]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW);
  Teuchos::RCP<Operators::Operator> pc_block11 = sub_pks_[1]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW);


  auto tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->PushBack(CreateTVSwithOneLeaf(op0->DomainMap()));
  tvs->PushBack(CreateTVSwithOneLeaf(op1->DomainMap()));

  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_matrix_->set_operator_block(0, 0, op0);
  op_tree_matrix_->set_operator_block(1, 1, op1);
  
  op_tree_pc_ = Teuchos::rcp(new Operators::TreeOperator(tvs));
  op_tree_pc_->set_operator_block(0, 0, pc_block00);
  op_tree_pc_->set_operator_block(1, 1, pc_block11);

  op_tree_rhs_ = Teuchos::rcp(new TreeVector());
  op_tree_rhs_->PushBack(CreateTVwithOneLeaf(op0->rhs()));
  op_tree_rhs_->PushBack(CreateTVwithOneLeaf(op1->rhs()));

  
  // full preconditioner
  if (precon_type_ == PRECON_FULL) {

    // set up the operator dWC_dT_block_
    Teuchos::ParameterList divq_plist(glist_->sublist("PKs").
                                      sublist(pks[0]).
                                      sublist("operators").
                                      sublist("diffusion operator").
                                      sublist("preconditioner"));

    // if (is_fv_) divq_plist.set("Newton correction", "true Jacobian");
    // else divq_plist.set("Newton correction", "approximate Jacobian");

    divq_plist.set("Newton correction", "true Jacobian");     
    divq_plist.set("exclude primary terms", true);

    Operators::PDE_DiffusionFactory opfactory;
    
    //ddivq_dT_ = opfactory.CreateWithGravity(divq_plist, mesh_);
    //dWC_dT_block_ = ddivq_dT_->global_operator();

    //ddivq_dT_->SetBCs(sub_pks_[0]->op_bc(), sub_pks_[1]->op_bc());
    
    if (dWC_dT_block_ == Teuchos::null) {
      dWC_dT_ =
        Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
      dWC_dT_block_ = dWC_dT_->global_operator();
    } else {
      dWC_dT_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, dWC_dT_block_));
    }


    // set up the operator dE_dp_block
    Teuchos::ParameterList ddivKgT_dp_plist(glist_->sublist("PKs").
                                      sublist(pks[1]).
                                      sublist("operators").
                                      sublist("diffusion operator").
                                      sublist("preconditioner"));

    // if (is_fv_) ddivKgT_dp_plist.set("Newton correction", "true Jacobian");
    // else ddivKgT_dp_plist.set("Newton correction", "approximate Jacobian");
    ddivKgT_dp_plist.set("Newton correction", "true Jacobian");
    ddivKgT_dp_plist.set("exclude primary terms", true);


    // ddivKgT_dp_ = opfactory.Create(ddivKgT_dp_plist, mesh_);
    // dE_dp_block_ = ddivKgT_dp_->global_operator();
      
    // ddivKgT_dp_->SetBCs(sub_pks_[1]->op_bc(), sub_pks_[0]->op_bc());


    if (dE_dp_block_ == Teuchos::null) {
      dE_dp_ =
        Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
      dE_dp_block_ = dE_dp_->global_operator();
    } else {
      dE_dp_ = Teuchos::rcp(
                            new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, dE_dp_block_));
    }
      
    Teuchos::ParameterList divhq_dp_plist(glist_->sublist("PKs").
                                          sublist(pks[1]).
                                          sublist("operators").
                                          sublist("diffusion operator").
                                          sublist("preconditioner"));

    // if (is_fv_) divhq_dp_plist.set("Newton correction", "true Jacobian");
    // else divhq_dp_plist.set("Newton correction", "approximate Jacobian");
    divhq_dp_plist.set("Newton correction", "true Jacobian");
    divhq_dp_plist.set("exclude primary terms", true);

    //Operators::PDE_DiffusionFactory opfactory;
    // if (dE_dp_block_ == Teuchos::null) {
    //   ddivhq_dp_ = opfactory.CreateWithGravity(divhq_dp_plist, mesh_);
    //   dE_dp_block_ = ddivhq_dp_->global_operator();
    // } else {
    //   ddivhq_dp_ = opfactory.CreateWithGravity(divhq_dp_plist, dE_dp_block_);
    // }

    /*--------------------*/
    // derivative with respect to temperature
    Teuchos::ParameterList divhq_dT_plist(glist_->sublist("PKs").
                                          sublist(pks[1]).
                                          sublist("operators").
                                          sublist("diffusion operator").
                                          sublist("preconditioner"));

    

    // if (is_fv_) divhq_dT_plist.set("Newton correction", "true Jacobian");
    // else divhq_dT_plist.set("Newton correction", "approximate Jacobian");
    // divhq_dT_plist.set("Newton correction", "true Jacobian");
    // divhq_dT_plist.set("exclude primary terms", true);
    
    // ddivhq_dT_ = opfactory.CreateWithGravity(divhq_dT_plist, pc_block11);

    op_tree_pc_->set_operator_block(0, 1, dWC_dT_block_);      
    op_tree_pc_->set_operator_block(1, 0, dE_dp_block_);
  }
  else if(precon_type_ == PRECON_FINITEDIFF){

    Teuchos::ParameterList divq_plist(glist_->sublist("PKs").
                                      sublist(pks[0]).
                                      sublist("operators").
                                      sublist("diffusion operator").
                                      sublist("preconditioner"));

    // if (is_fv_) divq_plist.set("Newton correction", "true Jacobian");
    // else divq_plist.set("Newton correction", "approximate Jacobian");

    divq_plist.set("Newton correction", "true Jacobian");     
    divq_plist.set("exclude primary terms", true);

    Operators::PDE_DiffusionFactory opfactory;
    
    ddivq_dT_ = opfactory.CreateWithGravity(divq_plist, mesh_);
    dWC_dT_block_ = ddivq_dT_->global_operator();
    //ddivq_dT_->SetBCs(sub_pks_[0]->op_bc(), sub_pks_[1]->op_bc());

    //dWC_dT_block_->PrintDiagnostics();    
    // dWC_dT_mat->PutScalar(0.0);
    // dWC_dT_mat->FillComplete();
    dWC_dT_block_->SymbolicAssembleMatrix();
    dWC_dT_block_->AssembleMatrix();
    Teuchos::RCP<Epetra_CrsMatrix> dWC_dT_mat = dWC_dT_block_->A();
    dWC_dT_mat->PutScalar(0.0);
        
    op_tree_pc_->set_operator_block(0, 1, dWC_dT_block_);      
    
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

  //exit(0);
  
}


/* *******************************************************************
* Performs one timestep.
******************************************************************* */
bool
FlowEnergy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
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
FlowEnergy_PK::FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f)
{

  double P0 = 2.10000000000000000e+07;
  double T0 = 700.0;

  //*vo_->os() << "Residual" << std::endl;
  std::vector<std::string> vnames;
  vnames.push_back("p");
  vnames.push_back("T");

  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(u_new->SubVector(0)->Data().ptr());
  vecs.push_back(u_new->SubVector(1)->Data().ptr());

  db_->WriteVectors(vnames, vecs, false);
  
  // flow
  auto u_old0 = u_old->SubVector(0);
  auto u_new0 = u_new->SubVector(0);

  double p_max, p_min;
  u_new0->Data()->MaxValue(&p_max);
  u_new0->Data()->MinValue(&p_min);
  //  std::cout<<"p_max="<<p_max<<" p_min="<<p_min<<"\n";

// energy
  auto u_old1 = u_old->SubVector(1);
  auto u_new1 = u_new->SubVector(1);
  auto f1 = f->SubVector(1);

  double t_max, t_min;
  u_new1->Data()->MaxValue(&t_max);
  u_new1->Data()->MinValue(&t_min);
  //u_new1->Data()->Print(*vo_->os());
  //std::cout<<"t_max="<<t_max<<" t_min="<<t_min<<"\n";
  
  
  auto f0 = f->SubVector(0);
  sub_pks_[0]->FunctionalResidual(t_old, t_new, u_old0, u_new0, f0);

  double f_norm;

  f0->Norm2(&f_norm);
  f0->Scale(1. / P0);
  f0->Norm2(&f_norm);

  // update molar flux
  Key key = Keys::getKey(domain_, "molar_flow_rate");
  auto mol_flowrate = S_->GetPtrW<CV_t>(key, Tags::DEFAULT, "");
  auto op0 = sub_pks_[0]->my_pde(Operators::PDEType::PDE_DIFFUSION);
  op0->UpdateFlux(u_new0->Data().ptr(), mol_flowrate.ptr());

  if (Keys::getVarName(sub_pks_[0]->name()) == "darcy") {
    double molar_mass = S_->Get<double>("const_fluid_molar_mass");
    mol_flowrate->Scale(1.0 / molar_mass);
  }

    
  sub_pks_[1]->FunctionalResidual(t_old, t_new, u_old1, u_new1, f1);


  f1->Norm2(&f_norm);
  f1->Scale(1./T0);
  f1->Norm2(&f_norm);

  
}


/* *******************************************************************
* Preconditioner update
******************************************************************* */
void
FlowEnergy_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt)
{
  PK_MPCStrong<PK_BDF>::UpdatePreconditioner(t, up, dt);

  Teuchos::RCP<Epetra_CrsMatrix> dWC_dT_mat = dWC_dT_block_->A();
  
  int ncells = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  if (precon_type_ == PRECON_FINITEDIFF) {

    WriteStateStatistics(*S_, *vo_);
    
    auto temp = *S_->GetPtrW<CompositeVector>(temperature_key_, Tags::DEFAULT, "")->ViewComponent("cell");
    auto pres = *S_->GetPtrW<CompositeVector>(pressure_key_, Tags::DEFAULT, "")->ViewComponent("cell");
    auto up_temp = *up->SubVector(1)->Data()->ViewComponent("cell");
    auto up_pres = *up->SubVector(0)->Data()->ViewComponent("cell");

    for (int n=0; n<ncells; n++){
      temp[0][n] = up_temp[0][n];
      pres[0][n] = up_pres[0][n];
    }
    sub_pks_[0]->ChangedSolution();
    sub_pks_[1]->ChangedSolution();        
    WriteStateStatistics(*S_, *vo_);

    auto u0 = Teuchos::rcp(new TreeVector(*up));
    u0->SubVector(0)->SetData(S_->GetPtrW<CompositeVector>(pressure_key_, Tags::DEFAULT, ""));
    u0->SubVector(1)->SetData(S_->GetPtrW<CompositeVector>(temperature_key_, Tags::DEFAULT, ""));
    
    auto u1 = Teuchos::rcp(new TreeVector(*u0));
    auto f0 = Teuchos::rcp(new TreeVector(*u0));
    auto f1 = Teuchos::rcp(new TreeVector(*u0));

  
    double tmax, factor, eps(1e-6), t_new, t_old;
    up->SubVector(1)->NormInf(&tmax);                                                                                                                                                                        
    factor = eps * tmax;
    t_old = t;
    t_new = t_old + dt;
  
    FunctionalResidual(t_old, t_new, up, up, f0);
  
    for (int n=0; n<ncells; n++){
      sub_pks_[0]->ChangedSolution();
      sub_pks_[1]->ChangedSolution();
      
      auto u1_tv = u1->SubVector(1);
      auto u1_cv = u1_tv->Data();
      
      auto& u1_v = *u1_cv->ViewComponent("cell");
      u1_v[0][n] += factor;
      for (int i=0; i<ncells; i++){
        std::cout<<"temp "<<u1_v[0][i]<<"\n";
      }
        
      FunctionalResidual(t_old, t_new, u0, u1, f1);
      u1_v[0][n] -= factor;

      auto& f0_v = *f0->SubVector(0)->Data()->ViewComponent("cell");
      auto& f1_v = *f1->SubVector(0)->Data()->ViewComponent("cell");


      for (int i=0; i<ncells; i++){
        std::cout<<"res "<<f1_v[0][i]<<" "<< f0_v[0][i]<<"\n";
      }
      
      int row_min(n), row_max(n);
      if (row_min>0) row_min--;
      if (row_max<ncells-1) row_max++;
      for (int r=row_min; r<=row_max; r++) {
        double val = (f1_v[0][r] - f0_v[0][r]) / factor;
        std::cout<< n<<" "<<r<<" "<<val<<"\n";
        dWC_dT_mat->InsertMyValues(r, 1, &val, &n);
      }      
    }
  }
  
  
  // if (include_pt_coupling_) {
  //   std::string passwd("");
  //   S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd, pressure_key_, Tags::DEFAULT);
  //   const auto& dEdP =
  //     S_->GetDerivative<CV_t>(energy_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT);

  //   S_->GetEvaluator(ws_key_).UpdateDerivative(*S_, passwd, temperature_key_, Tags::DEFAULT);
  //   const auto& dws_dT =
  //     S_->GetDerivative<CV_t>(ws_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);

  //   if (dt > 0.0) {
  //     op10_acc_->AddAccumulationDelta(*up->SubVector(0)->Data(), dEdP, dEdP, dt, "cell");
  //     op01_acc_->AddAccumulationDelta(*up->SubVector(1)->Data(), dws_dT, dws_dT, dt, "cell");
  //   }
  // }

  op_tree_pc_->AssembleMatrix();
  op_tree_pc_->InitializeInverse();
  op_tree_pc_->ComputeInverse();

  *vo_->os() << "UpdatePreconditioner." << std::endl;

  std::stringstream filename_s2;
  filename_s2 << "assembled_matrix" << 0 << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *dWC_dT_mat);

  
  std::vector<std::string> vnames;
  vnames.push_back("p");
  vnames.push_back("T");

  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(up->SubVector(0)->Data().ptr());
  vecs.push_back(up->SubVector(1)->Data().ptr());

  db_->WriteVectors(vnames, vecs, false);

  exit(0);
}


/* *******************************************************************
* Selection of default or full preconditioner
******************************************************************* */
int
FlowEnergy_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  int ierr;
  if (include_pt_coupling_) {
    Y->PutScalar(0.0);
    ierr = op_tree_pc_->ApplyInverse(*X, *Y);    
  } else {
    ierr = PK_MPCStrong<PK_BDF>::ApplyPreconditioner(X, Y);
  }

  *vo_->os() << "ApplyPreconditioner." << std::endl;
  std::vector<std::string> vnames;
  vnames.push_back("res p");
  vnames.push_back("PC*resp");
  vnames.push_back("res T");
  vnames.push_back("PC*resT");

  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(X->SubVector(0)->Data().ptr());
  vecs.push_back(Y->SubVector(0)->Data().ptr());
  vecs.push_back(X->SubVector(1)->Data().ptr());
  vecs.push_back(Y->SubVector(1)->Data().ptr());

  db_->WriteVectors(vnames, vecs, false);

  return ierr;
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
FlowEnergy_PK::ModifyCorrection(double h,
                                Teuchos::RCP<const TreeVector> res,
                                Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<TreeVector> du){

  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified =
    PK_MPCStrong<PK_BDF>::ModifyCorrection(h, res, u, du);

  Teuchos::RCP<const TreeVector> flow_u = u->SubVector(0);
  Teuchos::RCP<const TreeVector> flow_res = res->SubVector(0);
  Teuchos::RCP<TreeVector> flow_du = du->SubVector(0);

  Teuchos::RCP<const TreeVector> ener_u = u->SubVector(1);
  Teuchos::RCP<const TreeVector> ener_res = res->SubVector(1);
  Teuchos::RCP<TreeVector> ener_du = du->SubVector(1);

  double P0 = 2.10000000000000000e+07;
  double T0 = 700.0;

  flow_du -> Scale(P0);
  ener_du -> Scale(T0);
  

  return modified;

}

// -----------------------------------------------------------------------------
// L-scheme stability constant is updated by PKs.
// -----------------------------------------------------------------------------
std::vector<Key>
FlowEnergy_PK::SetupLSchemeKey(Teuchos::ParameterList& plist)
{
  auto tmp = sub_pks_[0]->SetupLSchemeKey(plist);
  L_scheme_keys_.insert(L_scheme_keys_.end(), tmp.begin(), tmp.end());
  return L_scheme_keys_;
}

} // namespace Amanzi
