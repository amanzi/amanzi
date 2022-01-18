/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  Process kernel for coupling Flow PK with Energy PK.
*/

#include "Energy_PK.hh"
#include "Flow_PK.hh"
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
                             const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // we will use a few parameter lists
  auto pk_list = Teuchos::sublist(glist, "PKs", true);
  my_list_ = Teuchos::sublist(pk_list, pk_name, true);
  domain_ = my_list_->template get<std::string>("domain name", "domain");
  
  Teuchos::ParameterList vlist;
  vo_ = Teuchos::rcp(new VerboseObject("FlowEnergy-" + domain_, vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void FlowEnergy_PK::Setup()
{
  mesh_ = S_->GetMesh(domain_);

  Teuchos::ParameterList& elist = S_->FEList();

  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(my_list_, "physical models and assumptions");
  bool vapor_diff = physical_models->get<bool>("vapor diffusion");

  if (physical_models->isParameter("eos lookup table")) {
    eos_table_ = physical_models->get<std::string>("eos lookup table");
  }

  // keys
  particle_density_key_ = Keys::getKey(domain_, "particle_density");
  ie_rock_key_ = Keys::getKey(domain_, "internal_energy_rock");
  ie_gas_key_ = Keys::getKey(domain_, "internal_energy_gas");
  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");

  energy_key_ = Keys::getKey(domain_, "energy");
  prev_energy_key_ = Keys::getKey(domain_, "prev_energy");

  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mol_density_gas_key_ = Keys::getKey(domain_, "molar_density_gas");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");

  pressure_key_ = Keys::getKey(domain_, "pressure");
  sat_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  prev_sat_liquid_key_ = Keys::getKey(domain_, "prev_saturation_liquid");

  wc_key_ = Keys::getKey(domain_, "water_content");
  prev_wc_key_ = Keys::getKey(domain_, "prev_water_content");

  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid");

  // Require primary field for this PK, which is pressure
  // Fields for solids
  // -- rock
  if (!S_->HasData(particle_density_key_)) {
    S_->Require<CV_t, CVS_t>(particle_density_key_, Tags::DEFAULT, particle_density_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(particle_density_key_);
  }

  if (!S_->HasData(ie_rock_key_) && !elist.isSublist(ie_rock_key_)) {
    Teuchos::Array<std::string> regions({ "All" });
    elist.sublist(ie_rock_key_)
         .set<std::string>("evaluator type", "iem")
         .set<std::string>("internal energy key", ie_rock_key_);
    elist.sublist(ie_rock_key_).sublist("IEM parameters").sublist("Material 1")
         .set<Teuchos::Array<std::string> >("regions", regions).sublist("IEM parameters")
         .set<std::string>("iem type", "linear")
         .set<double>("heat capacity", 620.0);
  }

  // Fields for gas
  // -- internal energy
  if (!S_->HasData(ie_gas_key_) && !elist.isSublist(ie_gas_key_)) {
    elist.sublist(ie_gas_key_)
         .set<std::string>("evaluator type", "iem water vapor")
         .set<std::string>("internal energy key", ie_gas_key_);
  }

  // -- molar density
  if (!S_->HasData(mol_density_gas_key_) && !elist.isSublist(mol_density_gas_key_)) {
    elist.sublist(mol_density_gas_key_)
         .set<std::string>("evaluator type", "eos")
         .set<std::string>("eos basis", "molar")
         .set<std::string>("molar density key", mol_density_gas_key_);
    elist.sublist(mol_density_gas_key_).sublist("EOS parameters")
         .set<std::string>("eos type", "vapor in gas");
    elist.sublist(mol_density_gas_key_).sublist("EOS parameters")
         .sublist("gas EOS parameters")
         .set<std::string>("eos type", "ideal gas")
         .set<double>("molar mass of gas", 28.9647e-03);  // dry air
  }

  // -- molar fraction FIXME (it is not used by all models)
  if (!S_->HasData("molar_fraction_gas") && !elist.isSublist("molar_fraction_gas")) {
    elist.sublist("molar_fraction_gas")
         .set<std::string>("evaluator type", "molar fraction gas")
         .set<std::string>("molar fraction key", "molar_fraction_gas");
    elist.sublist("molar_fraction_gas")
         .sublist("vapor pressure model parameters")
         .set<std::string>("eos type", "water vapor over water/ice");
  }

  // Fields for liquid
  // -- internal energy
  S_->Require<CV_t, CVS_t>(ie_liquid_key_, Tags::DEFAULT, ie_liquid_key_)
    .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(ie_liquid_key_);


  // -- molar and mass density
  S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
    .SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(mol_density_liquid_key_);

  S_->RequireEvaluator(mass_density_liquid_key_);
  if (S_->GetEvaluator(mass_density_liquid_key_).IsDifferentiableWRT(*S_, pressure_key_, Tags::DEFAULT)) {
    S_->RequireDerivative<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT,
                                      pressure_key_, Tags::DEFAULT, mass_density_liquid_key_);
  }

  // -- viscosity
  S_->Require<CV_t, CVS_t>(viscosity_liquid_key_, Tags::DEFAULT, viscosity_liquid_key_)
    .SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(viscosity_liquid_key_);

  // inform other PKs about strong coupling
  // -- flow
  auto pks = glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string> >("PKs order").toVector();
  std::string model = (vapor_diff) ? "two-phase" : "one-phase";
  Teuchos::ParameterList& flow = glist_->sublist("PKs").sublist(pks[0])
                                        .sublist("physical models and assumptions");
  flow.set<bool>("vapor diffusion", vapor_diff);

  // -- energy
  Teuchos::ParameterList& energy = glist_->sublist("PKs").sublist(pks[1])
                                          .sublist("physical models and assumptions");
  energy.set<bool>("vapor diffusion", vapor_diff);

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup();

  // copies of fields (must be called after PKs)
  if (S_->HasData(prev_wc_key_)) {
    S_->Require<CV_t, CVS_t>(prev_wc_key_, Tags::COPY, prev_wc_key_);
    S_->GetRecordW(prev_wc_key_, Tags::COPY, prev_wc_key_).set_initialized();
  }
}


/* ******************************************************************* 
* Initialization of copies requires fileds to exists
******************************************************************* */
void FlowEnergy_PK::Initialize()
{
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
  op_tree_pc_->set_operator_block(0, 0, sub_pks_[0]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW));
  op_tree_pc_->set_operator_block(1, 1, sub_pks_[1]->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW));

  op_tree_rhs_ = Teuchos::rcp(new TreeVector());
  op_tree_rhs_->PushBack(CreateTVwithOneLeaf(op0->rhs()));
  op_tree_rhs_->PushBack(CreateTVwithOneLeaf(op1->rhs()));

  // output of initialization statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
               << "matrix:" << std::endl
               << op_tree_matrix_->PrintDiagnostics() << std::endl
               << "preconditioner:" << std::endl
               << op_tree_pc_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt() 
               << vo_->reset() << std::endl << std::endl;
  }

  AMANZI_ASSERT(sub_pks_.size() == 2);
}
  

/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowEnergy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // flow
  // -- swap saturations (current and previous)
  S_->GetEvaluator(sat_liquid_key_).Update(*S_, "flow");
  const auto& sat = S_->Get<CV_t>(sat_liquid_key_);
  auto& sat_prev = S_->GetW<CV_t>(prev_sat_liquid_key_, Tags::DEFAULT, "flow");

  CompositeVector sat_prev_copy(sat_prev);
  sat_prev = sat;

  // -- swap water_contents (current and previous)
  if (S_->HasData(wc_key_)) {
    S_->Copy(prev_wc_key_, Tags::COPY, Tags::DEFAULT);

    S_->GetEvaluator(wc_key_).Update(*S_, "flow");
    const auto& wc = S_->Get<CV_t>(wc_key_);
    auto& wc_prev = S_->GetW<CV_t>(prev_wc_key_, "flow");
    wc_prev = wc;
  }

  // energy
  // -- swap conserved energies (current and previous)
  S_->GetEvaluator(energy_key_).Update(*S_, "thermal");
  const auto& e = S_->Get<CV_t>(energy_key_);
  auto& e_prev = S_->GetW<CV_t>(prev_energy_key_, "thermal");

  CompositeVector e_prev_copy(e_prev);
  e_prev = e;
 
  // try time step
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    // recover conserved quantaties at the beginning of time step
    S_->GetW<CV_t>(prev_sat_liquid_key_, "flow") = sat_prev_copy;
    S_->GetW<CV_t>(prev_energy_key_, "thermal") = e_prev_copy;
    if (S_->HasData(wc_key_)) {
      S_->Copy(prev_wc_key_, Tags::DEFAULT, Tags::COPY);
    }

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed. Restored " << prev_sat_liquid_key_ 
               << ", " << prev_wc_key_ << ", " << prev_energy_key_ << std::endl;
  }

  return fail;
}

}  // namespace Amanzi

