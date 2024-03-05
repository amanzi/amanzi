/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Iterative coupling of mechanics and flow PKs via undrained split.
*/

#include <string>

#include "PK_Physical.hh"
#include "PK_Utils.hh"
#include "StateArchive.hh"

#include "MechanicsFlow_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Constructor
****************************************************************** */
MechanicsFlow_PK::MechanicsFlow_PK(Teuchos::ParameterList& pk_tree,
                                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCSequential(pk_tree, glist, S, soln), glist_(glist), thermal_flow_(false)
{
  domain_ = my_list_->template get<std::string>("domain name", "domain");

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = my_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("MechanicsFlow", vlist));
}


/* ******************************************************************
* Setup of PK
****************************************************************** */
void
MechanicsFlow_PK::Setup()
{
  std::string passwd("");

  displacement_key_ = Keys::getKey(domain_, "displacement"); // primaries
  pressure_key_ = Keys::getKey(domain_, "pressure");

  hydrostatic_stress_key_ = Keys::getKey(domain_, "hydrostatic_stress");
  vol_strain_key_ = Keys::getKey(domain_, "volumetric_strain");
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  water_storage_key_ = Keys::getKey(domain_, "water_storage");

  porosity_key_ = Keys::getKey(domain_, "porosity");
  undrained_split_coef_key_ = Keys::getKey(domain_, "undrained_split_coef");

  S_->Require<CV_t, CVS_t>(hydrostatic_stress_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CV_t, CVS_t>(vol_strain_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CV_t, CVS_t>(undrained_split_coef_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  for (int i = 0; i < 2; ++i) {
    glist_->sublist("PKs")
      .sublist(pks[i])
      .sublist("physical models and assumptions")
      .set<bool>("biot scheme: undrained split", true)
      .set<bool>("biot scheme: fixed stress split", false)
      .set<bool>("thermoelasticity", thermal_flow_);
  }

  PK_MPCSequential::Setup();
}


/* ******************************************************************
* Extend default initialization algorithms
****************************************************************** */
void
MechanicsFlow_PK::Initialize()
{
  InitializeCVField(S_, *vo_, undrained_split_coef_key_, Tags::DEFAULT, "", 0.0);
  PK_MPCSequential::Initialize();
}


/* ******************************************************************
* Extended treatment of time step in transport PK.
****************************************************************** */
bool
MechanicsFlow_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  std::vector<std::string> fields(
    { pressure_key_, saturation_liquid_key_, displacement_key_, vol_strain_key_ });
  if (S_->HasRecord(water_storage_key_)) fields.push_back(water_storage_key_);

  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

  // compute stability coefficient
  if (Keys::getVarName(sub_pks_[1]->name()) == "darcy") {
    Exceptions::amanzi_throw("Saturated flow model is under development. Use Richards' model.");
    EvaluateForDarcy_();
  } else {
    EvaluateForRichards_();
  }

  bool fail = PK_MPCSequential::AdvanceStep(t_old, t_new, reinit);
  if (fail) archive.Restore("");

  return fail;
}


/* ******************************************************************
* Stability coefficient for undrained split algorithm
****************************************************************** */
void
MechanicsFlow_PK::EvaluateForDarcy_()
{
  Key specific_storage_key = Keys::getKey(domain_, "specific_storage");
  const auto& ss_c = *S_->Get<CompositeVector>(specific_storage_key).ViewComponent("cell");

  int d = S_->GetMesh(domain_)->getSpaceDimension();
  auto gravity = S_->Get<AmanziGeometry::Point>("gravity");
  double g = fabs(gravity[d - 1]);
  double rho = S_->Get<double>("const_fluid_density");

  auto tmp = const_cast<Evaluator*>(&S_->GetEvaluator(porosity_key_));
  auto eval = dynamic_cast<Flow::PorosityEvaluator*>(tmp);

  auto& stability_c =
    *S_->GetW<CompositeVector>(undrained_split_coef_key_, "").ViewComponent("cell");
  int ncells = stability_c.MyLength();

  double stbmax(0.0);
  for (int c = 0; c != ncells; ++c) {
    double b = eval->getBiotCoefficient(c);
    stability_c[0][c] = rho * g * b * b / ss_c[0][c];
    stbmax = std::max(stbmax, stability_c[0][c]);
  }

  if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "max of stabilization: " << stbmax << std::endl;
  }
}

void
MechanicsFlow_PK::EvaluateForRichards_()
{
  Key mol_density_key = Keys::getKey(domain_, "molar_density_liquid");
  const auto& eta_c = *S_->Get<CompositeVector>(mol_density_key).ViewComponent("cell");

  S_->GetEvaluator(water_storage_key_).UpdateDerivative(*S_, "mpc", pressure_key_, Tags::DEFAULT);
  const auto& dws_dp = *S_->GetDerivative<CompositeVector>(
                            water_storage_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT)
                          .ViewComponent("cell");

  auto tmp = const_cast<Evaluator*>(&S_->GetEvaluator(porosity_key_));
  auto eval = dynamic_cast<Flow::PorosityEvaluator*>(tmp);

  auto& stability_c =
    *S_->GetW<CompositeVector>(undrained_split_coef_key_, "").ViewComponent("cell");
  int ncells = stability_c.MyLength();

  double stbmax(0.0);
  for (int c = 0; c != ncells; ++c) {
    double b = eval->getBiotCoefficient(c);
    stability_c[0][c] = eta_c[0][c] * b * b / std::fabs(dws_dp[0][c]);
    stbmax = std::max(stbmax, stability_c[0][c]);
  }

  if (vo_->getVerbLevel() > Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "max of stabilization: " << stbmax << std::endl;
  }
}

} // namespace Amanzi
