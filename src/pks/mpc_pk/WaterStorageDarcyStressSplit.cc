/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MPC PK

  Field evaluator for water storage with a correction for fixed
  stress split.
*/

#include "CommonDefs.hh"

#include "WaterStorageDarcyStressSplit.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor.
****************************************************************** */
WaterStorageDarcyStressSplit::WaterStorageDarcyStressSplit(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("water storage key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  pressure_key_ = Keys::getKey(domain, "pressure");
  specific_storage_key_ = Keys::getKey(domain, "specific_storage");

  young_modulus_key_ = Keys::getKey(domain, "young_modulus");
  poisson_ratio_key_ = Keys::getKey(domain, "poisson_ratio");
  biot_key_ = Keys::getKey(domain, "biot_coefficient");
  strain_key_ = Keys::getKey(domain, "volumetric_strain");

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(strain_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
WaterStorageDarcyStressSplit::WaterStorageDarcyStressSplit(
  const WaterStorageDarcyStressSplit& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other){};


Teuchos::RCP<Evaluator>
WaterStorageDarcyStressSplit::Clone() const
{
  return Teuchos::rcp(new WaterStorageDarcyStressSplit(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageDarcyStressSplit::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& results)
{
  const auto& p = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& ss = *S.Get<CompositeVector>(specific_storage_key_).ViewComponent("cell");

  double rho = S.Get<double>("const_fluid_density");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  const auto& E = *S.Get<CompositeVector>(young_modulus_key_).ViewComponent("cell");
  const auto& nu = *S.Get<CompositeVector>(poisson_ratio_key_).ViewComponent("cell");
  const auto& b = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");
  const auto& e = *S.Get<CompositeVector>(strain_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    double stability = FixedStressStabilityDarcy(E[0][c], nu[0][c], b[0][c]);
    result_v[0][c] = p[0][c] * ss[0][c] / g + b[0][c] * e[0][c] + stability * rho * p[0][c];
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageDarcyStressSplit::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  const auto& ss = *S.Get<CompositeVector>(specific_storage_key_).ViewComponent("cell");
  double rho = S.Get<double>("const_fluid_density");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  const auto& E = *S.Get<CompositeVector>(young_modulus_key_).ViewComponent("cell");
  const auto& nu = *S.Get<CompositeVector>(poisson_ratio_key_).ViewComponent("cell");
  const auto& b = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      double stability = FixedStressStabilityDarcy(E[0][c], nu[0][c], b[0][c]);
      result_v[0][c] = ss[0][c] / g + stability * rho;
    }
  } else if (wrt_key == strain_key_) {
    for (int c = 0; c != ncells; ++c) result_v[0][c] = b[0][c];
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Amanzi
