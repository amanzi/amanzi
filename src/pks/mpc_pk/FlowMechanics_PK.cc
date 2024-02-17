/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Weak coupling of mechanics and flow PKs.
*/

#include <string>

#include "PK_BDF.hh"
#include "PorosityEvaluator.hh"
#include "StateArchive.hh"
#include "Transport_PK.hh"

#include "FlowMechanics_PK.hh"
#include "WaterStorageStressSplit.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

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
  displacement_key_ = Keys::getKey(domain_, "displacement"); // primary
  hydrostatic_stress_key_ = Keys::getKey(domain_, "hydrostatic_stress");
  vol_strain_key_ = Keys::getKey(domain_, "volumetric_strain");

  pressure_key_ = Keys::getKey(domain_, "pressure"); // primary
  porosity_key_ = Keys::getKey(domain_, "porosity");
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  water_storage_key_ = Keys::getKey(domain_, "water_storage");

  thermal_flow_ = (Keys::getVarName(sub_pks_[0]->name()) == "thermal flow");

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

  // flow
  // -- we redefine ws to get the fixed stress split scheme
  S_->Require<CV_t, CVS_t>(water_storage_key_, Tags::DEFAULT, water_storage_key_)
    .SetMesh(mesh)
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::ParameterList elist(water_storage_key_);
  elist.set<std::string>("water storage key", water_storage_key_)
    .set<std::string>("tag", "")
    .set<std::string>("pressure key", pressure_key_)
    .set<std::string>("saturation key", saturation_liquid_key_)
    .set<std::string>("porosity key", porosity_key_);

  S_->RequireDerivative<CV_t, CVS_t>(
      water_storage_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT, water_storage_key_)
    .SetGhosted();

  auto eval = Teuchos::rcp(new WaterStorageStressSplit(elist));
  S_->SetEvaluator(water_storage_key_, Tags::DEFAULT, eval);

  PK_MPCSequential::Setup();
}


/* ******************************************************************
* Overloading time-stepping to capture failure
****************************************************************** */
bool
FlowMechanics_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  std::vector<std::string> fields({ pressure_key_,
                                    saturation_liquid_key_,
                                    water_storage_key_,
                                    displacement_key_,
                                    vol_strain_key_ });
  StateArchive archive(S_, vo_);
  archive.Add(fields, Tags::DEFAULT);

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
  Key prev = Keys::getKey(domain_, "prev_water_storage");
  auto& ws_c = *S_->GetW<CV_t>(prev, Tags::DEFAULT, "").ViewComponent("cell");

  // access to pressures, depends on PK
  Teuchos::RCP<const Epetra_MultiVector> u0_c, u1_c;
  if (thermal_flow_) {
    u0_c = u_old->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell");
    u1_c = u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell");
  } else {
    u0_c = u_old->SubVector(0)->Data()->ViewComponent("cell");
    u1_c = u_new->SubVector(0)->Data()->ViewComponent("cell");
  }

  const auto& E = *S_->Get<CompositeVector>("young_modulus").ViewComponent("cell");
  const auto& nu = *S_->Get<CompositeVector>("poisson_ratio").ViewComponent("cell");
  const auto& n_l = *S_->Get<CompositeVector>("molar_density_liquid").ViewComponent("cell");

  auto tmp = const_cast<Evaluator*>(&S_->GetEvaluator(porosity_key_));
  auto eval = dynamic_cast<Flow::PorosityEvaluator*>(tmp);

  int ncells = ws_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    double b = eval->getBiotCoefficient(c);
    double stability = FixedStressStability(E[0][c], nu[0][c], b);
    ws_c[0][c] += stability * n_l[0][c] * ((*u1_c)[0][c] - (*u0_c)[0][c]);
  }
}

} // namespace Amanzi