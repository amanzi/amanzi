/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MPC PK

  Field evaluator for water storage with a correction for volumetric strain.
*/

#include "CommonDefs.hh"

#include "WaterStorageDarcyPoroelasticity.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor.
****************************************************************** */
WaterStorageDarcyPoroelasticity::WaterStorageDarcyPoroelasticity(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("water storage key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  pressure_key_ = Keys::getKey(domain, "pressure");
  specific_storage_key_ = Keys::getKey(domain, "specific_storage");

  biot_key_ = Keys::getKey(domain, "biot_coefficient");
  strain_key_ = Keys::getKey(domain, "volumetric_strain");

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(strain_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
WaterStorageDarcyPoroelasticity::WaterStorageDarcyPoroelasticity(
  const WaterStorageDarcyPoroelasticity& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other) {};


Teuchos::RCP<Evaluator>
WaterStorageDarcyPoroelasticity::Clone() const
{
  return Teuchos::rcp(new WaterStorageDarcyPoroelasticity(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageDarcyPoroelasticity::Evaluate_(const State& S,
                                        const std::vector<CompositeVector*>& results)
{
  const auto& p = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& ss = *S.Get<CompositeVector>(specific_storage_key_).ViewComponent("cell");

  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  const auto& b = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");
  const auto& e = *S.Get<CompositeVector>(strain_key_).ViewComponent("cell");

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = p[0][c] * ss[0][c] / g + b[0][c] * e[0][c];
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageDarcyPoroelasticity::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  const auto& ss = *S.Get<CompositeVector>(specific_storage_key_).ViewComponent("cell");
  const auto& b = *S.Get<CompositeVector>(biot_key_).ViewComponent("cell");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = ss[0][c] / g;
    }
  } else if (wrt_key == strain_key_) {
    for (int c = 0; c != ncells; ++c) result_v[0][c] = b[0][c];
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Amanzi
