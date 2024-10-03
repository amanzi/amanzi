/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#include "errors.hh"

#include "PressureHeadEvaluator.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Simple constructor
****************************************************************** */
PressureHeadEvaluator::PressureHeadEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(
      std::make_pair(plist_.get<std::string>("my key", "pressure_head"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  diameter_key_ = plist_.get<std::string>("diameter key", Keys::getKey(domain, "diameter"));
  wetted_angle_key_ =
    plist_.get<std::string>("wetted angle key", Keys::getKey(domain, "wetted_angle"));
  wetted_area_key_ =
    plist_.get<std::string>("wetted area key", Keys::getKey(domain, "wetted_area"));

  celerity_ = plist_.get<double>("celerity", 2);

  dependencies_.insert(std::make_pair(wetted_area_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
PressureHeadEvaluator::Clone() const
{
  return Teuchos::rcp(new PressureHeadEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PressureHeadEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& WettedAngle_c = *S.Get<CompositeVector>(wetted_angle_key_).ViewComponent("cell");
  const auto& WettedArea_c = *S.Get<CompositeVector>(wetted_area_key_).ViewComponent("cell");
  const auto& PipeD_c = *S.Get<CompositeVector>(diameter_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    double pipeCrossSection = Pi * 0.25 * PipeD_c[0][c] * PipeD_c[0][c];
    if (WettedAngle_c[0][c] >= TwoPi) {
      result_c[0][c] =
        (celerity_ * celerity_ * (WettedArea_c[0][c] - pipeCrossSection)) / (g * pipeCrossSection);
    } else {
      result_c[0][c] = 0.0;
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
PressureHeadEvaluator::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& results)
{
  const auto& WettedAngle_c = *S.Get<CompositeVector>(wetted_angle_key_).ViewComponent("cell");
  const auto& PipeD_c = *S.Get<CompositeVector>(diameter_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  double g = norm(S.Get<AmanziGeometry::Point>("gravity"));

  int ncells = result_c.MyLength();
  if (wrt_key == wetted_angle_key_) {
    for (int c = 0; c != ncells; ++c) {
      double pipeCrossSection = Pi * 0.25 * PipeD_c[0][c] * PipeD_c[0][c];
      if (WettedAngle_c[0][c] >= TwoPi) {
        result_c[0][c] = (celerity_ * celerity_) / (g * pipeCrossSection);
        ;
      } else {
        result_c[0][c] = 0.0;
      }
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

} // namespace ShallowWater
} // namespace Amanzi
