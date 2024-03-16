/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Field evaluator for water storage which is the conserved quantity
  in the Darcy equation.
*/

#include "CommonDefs.hh"
#include "WaterStorageDarcy.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor.
****************************************************************** */
WaterStorageDarcy::WaterStorageDarcy(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), aperture_(false)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("water storage key"), Tags::DEFAULT));
  }
  pressure_key_ = plist_.get<std::string>("pressure key");
  specific_storage_key_ = plist_.get<std::string>("specific storage key");

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));

  if (plist_.isParameter("aperture key")) {
    aperture_ = true;
    aperture_key_ = plist_.get<std::string>("aperture key");
    dependencies_.insert(std::make_pair(aperture_key_, Tags::DEFAULT));
  }
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
WaterStorageDarcy::WaterStorageDarcy(const WaterStorageDarcy& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    aperture_(other.aperture_){};


Teuchos::RCP<Evaluator>
WaterStorageDarcy::Clone() const
{
  return Teuchos::rcp(new WaterStorageDarcy(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageDarcy::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& p = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& ss = *S.Get<CompositeVector>(specific_storage_key_).ViewComponent("cell");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  for (int c = 0; c != ncells; ++c) result_v[0][c] = p[0][c] * ss[0][c] / g;

  if (aperture_) {
    const auto& aperture = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) result_v[0][c] *= aperture[0][c];
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void
WaterStorageDarcy::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& results)
{
  const auto& p = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  const auto& ss = *S.Get<CompositeVector>(specific_storage_key_).ViewComponent("cell");
  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  auto& result_v = *results[0]->ViewComponent("cell");
  int ncells = results[0]->size("cell");

  if (wrt_key == pressure_key_) {
    for (int c = 0; c != ncells; ++c) result_v[0][c] = ss[0][c] / g;
  } else if (wrt_key == aperture_key_ && aperture_) {
    for (int c = 0; c != ncells; ++c) { result_v[0][c] = ss[0][c] * p[0][c] / g; }
  } else {
    AMANZI_ASSERT(0);
  }

  if (aperture_ && wrt_key != aperture_key_) {
    const auto& aperture = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");
    for (int c = 0; c != ncells; ++c) { result_v[0][c] *= aperture[0][c]; }
  }
}

} // namespace Flow
} // namespace Amanzi
