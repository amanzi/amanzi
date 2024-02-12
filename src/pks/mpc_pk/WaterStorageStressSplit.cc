/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  MPC PK

  Field evaluator for water storage with a correction for fixed
  stress split.
*/

#include "CommonDefs.hh"
#include "PorosityEvaluator.hh"

#include "WaterStorageStressSplit.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor.
****************************************************************** */
WaterStorageStressSplit::WaterStorageStressSplit(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("water storage key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  saturation_key_ = plist_.get<std::string>("saturation key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  mol_density_liquid_key_ = Keys::getKey(domain, "molar_density_liquid");
  pressure_key_ = Keys::getKey(domain, "pressure");

  young_modulus_key_ = Keys::getKey(domain, "young_modulus");
  poisson_ratio_key_ = Keys::getKey(domain, "poisson_ratio");

  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(saturation_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(mol_density_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
WaterStorageStressSplit::WaterStorageStressSplit(const WaterStorageStressSplit& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other){};


Teuchos::RCP<Evaluator>
WaterStorageStressSplit::Clone() const
{
  return Teuchos::rcp(new WaterStorageStressSplit(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageStressSplit::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& s_l = *S.Get<CompositeVector>(saturation_key_).ViewComponent("cell");
  const auto& n_l = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& p_l = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");

  const auto& E = *S.Get<CompositeVector>(young_modulus_key_).ViewComponent("cell");
  const auto& nu = *S.Get<CompositeVector>(poisson_ratio_key_).ViewComponent("cell");

  auto tmp = const_cast<Evaluator*>(&S.GetEvaluator(porosity_key_));
  auto eval = dynamic_cast<Flow::PorosityEvaluator*>(tmp);

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    double b = eval->getBiotCoefficient(c);
    double stability = FixedStressStability(E[0][c], nu[0][c], b);
    result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] + stability * n_l[0][c] * p_l[0][c];
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageStressSplit::EvaluatePartialDerivative_(const State& S,
                                                    const Key& wrt_key,
                                                    const Tag& wrt_tag,
                                                    const std::vector<CompositeVector*>& results)
{
  const auto& s_l = *S.Get<CompositeVector>(saturation_key_).ViewComponent("cell");
  const auto& n_l = *S.Get<CompositeVector>(mol_density_liquid_key_).ViewComponent("cell");
  const auto& phi = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");

  const auto& E = *S.Get<CompositeVector>(young_modulus_key_).ViewComponent("cell");
  const auto& nu = *S.Get<CompositeVector>(poisson_ratio_key_).ViewComponent("cell");

  auto tmp = const_cast<Evaluator*>(&S.GetEvaluator(porosity_key_));
  auto eval = dynamic_cast<Flow::PorosityEvaluator*>(tmp);

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = s_l[0][c] * n_l[0][c]; }
  } else if (wrt_key == saturation_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = phi[0][c] * n_l[0][c]; }
  } else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = phi[0][c] * s_l[0][c]; }
  } else if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      double b = eval->getBiotCoefficient(c);
      double stability = FixedStressStability(E[0][c], nu[0][c], b);
      result_v[0][c] = stability * n_l[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Amanzi
