/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Weak coupling of mechanics and flow PKs using the fixed-stress operator split.
*/

#include <string>

#include "PK_BDF.hh"
#include "StateArchive.hh"
#include "StateHelpers.hh"
#include "Transport_PK.hh"

#include "FlowMechanics_PK.hh"
#include "WaterStorageDarcyPoroelasticity.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

// Optimal stability coefficient due to Makilic, Wheeler.
// Convergence of iterative coupling for coupled flow and geomechanics.
// Comput Geosci 2013.
inline double
FixedStressStability(double E, double nu, double b)
{
  double mu = E / (2 * (1 + nu));
  return b * b / mu / 2;
}


/* ******************************************************************
* Constructor
****************************************************************** */
FlowMechanics_PK::FlowMechanics_PK(Teuchos::ParameterList& pk_tree,
                                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCSequential(pk_tree, glist, S, soln), glist_(glist), thermal_flow_(false)
{
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = my_list_->sublist("verbose object");

  domain_ = my_list_->template get<std::string>("domain name", "domain");
  vo_ = Teuchos::rcp(new VerboseObject("FlowMechanics", vlist));
}


/* ******************************************************************
* Setup of PK
****************************************************************** */
void
FlowMechanics_PK::Setup()
{
  std::string passwd("");
  pressure_key_ = Keys::getKey(domain_, "pressure");         // primary
  displacement_key_ = Keys::getKey(domain_, "displacement"); // primary

  hydrostatic_stress_key_ = Keys::getKey(domain_, "hydrostatic_stress");
  vol_strain_key_ = Keys::getKey(domain_, "volumetric_strain");
  biot_key_ = Keys::getKey(domain_, "biot_coefficient");

  porosity_key_ = Keys::getKey(domain_, "porosity");
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  water_storage_key_ = Keys::getKey(domain_, "water_storage");

  thermal_flow_ = (Keys::getVarName(sub_pks_[0]->name()) == "flow and energy");

  // mechanics
  auto mesh = S_->GetMesh(domain_);
  S_->Require<CV_t, CVS_t>(hydrostatic_stress_key_, Tags::DEFAULT, passwd)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CV_t, CVS_t>(vol_strain_key_, Tags::DEFAULT, passwd)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CV_t, CVS_t>(biot_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  for (int i = 0; i < 2; ++i) {
    glist_->sublist("PKs")
      .sublist(pks[i])
      .sublist("physical models and assumptions")
      .set<bool>("biot scheme: undrained split", false)
      .set<bool>("biot scheme: fixed stress split", true)
      .set<bool>("thermoelasticity", thermal_flow_);
  }

  // -- we re-define water_storage as a way to get the undrained split scheme
  if (Keys::getVarName(sub_pks_[0]->name()) == "darcy") {
    auto elist = RequireFieldForEvaluator(*S_, water_storage_key_);
    elist.set<std::string>("pressure key", pressure_key_);

    auto eval = Teuchos::rcp(new WaterStorageDarcyPoroelasticity(elist));
    S_->SetEvaluator(water_storage_key_, Tags::DEFAULT, eval);
  }

  S_->RequireDerivative<CV_t, CVS_t>(
      water_storage_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, water_storage_key_)
    .SetGhosted();

  PK_MPCSequential::Setup();
}


/* ******************************************************************
* Process additional parameter
****************************************************************** */
void
FlowMechanics_PK::Initialize()
{
  PK_MPCSequential::Initialize();

  if (my_list_->isParameter("initialize displacement")) {
    bool fail = sub_pks_[1]->AdvanceStep(0.0, 0.0, false);
    if (fail) Exceptions::amanzi_throw("Initialization of displacement has failed.");
    Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
      S_->GetEvaluatorPtr(displacement_key_, Tags::DEFAULT))
      ->SetChanged();

    S_->GetEvaluator(water_storage_key_).Update(*S_, "mpc");
    S_->GetW<CV_t>("prev_water_storage", Tags::DEFAULT, "") =
      S_->Get<CV_t>(water_storage_key_, Tags::DEFAULT);
  }
}


/* ******************************************************************
* Overloading time-stepping to capture failure
****************************************************************** */
bool
FlowMechanics_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  std::vector<std::string> fields(
    { pressure_key_, saturation_liquid_key_, displacement_key_, vol_strain_key_ });
  if (S_->HasRecord(water_storage_key_) ) fields.push_back(water_storage_key_);

  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

  // populate L-scheme stability constant
  if (L_scheme_) {
    Key key_stab = Keys::getKey(domain_, L_scheme_keys_[0] + "_stability");
    auto& stability_c = *S_->GetW<CV_t>(key_stab, "state").ViewComponent("cell");

    Key density_key = (Keys::getVarName(sub_pks_[0]->name()) == "darcy") ? "mass_density_liquid"
                                                                         : "molar_density_liquid";
    const auto& eta_c = *S_->Get<CompositeVector>(density_key).ViewComponent("cell");

    const auto& E = *S_->Get<CompositeVector>("young_modulus").ViewComponent("cell");
    const auto& nu = *S_->Get<CompositeVector>("poisson_ratio").ViewComponent("cell");
    const auto& b = *S_->Get<CompositeVector>("biot_coefficient").ViewComponent("cell");

    int ncells = b.MyLength();
    for (int c = 0; c != ncells; ++c) {
      double fss = FixedStressStability(E[0][c], nu[0][c], b[0][c]);
      stability_c[0][c] = fss * eta_c[0][c];
    }

    Key key_prev = Keys::getKey(domain_, L_scheme_keys_[0] + "_prev");
    *S_->GetW<CV_t>(key_prev, "state").ViewComponent("cell") = *S_->Get<CV_t>(pressure_key_).ViewComponent("cell");
  }

  bool fail = PK_MPCSequential::AdvanceStep(t_old, t_new, reinit);
  if (fail) archive.Restore("");

  return fail;
}


/* ******************************************************************
* Use error norms for each PK.
****************************************************************** */
double
FlowMechanics_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double err(0.0);
  for (int i = 0; i < 2; ++i) {
    auto pk = Teuchos::rcp_dynamic_cast<PK_BDF>(sub_pks_[i]);
    err = std::max(err, pk->ErrorNorm(u->SubVector(i), du->SubVector(i)));
  }
  return err;
}


/* ******************************************************************
* Update previous accumulation field
****************************************************************** */
void
FlowMechanics_PK::CommitSequentialStep(Teuchos::RCP<const TreeVector> u_old,
                                       Teuchos::RCP<const TreeVector> u_new)
{
  // access to pressures, depends on PK
  std::string name;
  Teuchos::RCP<const Epetra_MultiVector> u1_c;
  if (thermal_flow_) {
    auto mpc = Teuchos::rcp_dynamic_cast<PK_MPC<PK_BDF>>(sub_pks_[0]);
    name = Keys::getVarName((*mpc->begin())->name());
    u1_c = u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell");
  } else {
    name = Keys::getVarName(sub_pks_[0]->name());
    u1_c = u_new->SubVector(0)->Data()->ViewComponent("cell");
  }

  // swap current/previous solutions for L-scheme
  if (L_scheme_) {
    Key key_prev = Keys::getKey(domain_, L_scheme_keys_[0] + "_prev");
    *S_->GetW<CV_t>(key_prev, "state").ViewComponent("cell") = *u1_c;
  }
}

} // namespace Amanzi
