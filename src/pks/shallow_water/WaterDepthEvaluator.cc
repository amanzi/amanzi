/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

/*
  Pipe Flow PK

*/

#include "errors.hh"

#include "WaterDepthEvaluator.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Simple constructor
****************************************************************** */
WaterDepthEvaluator::WaterDepthEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(
      std::make_pair(plist_.get<std::string>("my key", "water_depth"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  diameter_key_ = plist_.get<std::string>("diameter key", Keys::getKey(domain, "diameter"));

  wetted_angle_key_ =
    plist_.get<std::string>("wetted angle key", Keys::getKey(domain, "wetted_angle"));

  primary_variable_key_ =
    plist_.get<std::string>("wetted area key", Keys::getKey(domain, "wetted_area"));

  dependencies_.insert(std::make_pair(wetted_angle_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(primary_variable_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
WaterDepthEvaluator::Clone() const
{
  return Teuchos::rcp(new WaterDepthEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
WaterDepthEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& WettedAngle_c = *S.Get<CompositeVector>(wetted_angle_key_).ViewComponent("cell");
  const auto& WettedArea_c = *S.Get<CompositeVector>(primary_variable_key_).ViewComponent("cell");
  const auto& PipeD_c = *S.Get<CompositeVector>(diameter_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    if (WettedAngle_c[0][c] >= (2.0 * M_PI)) {
      // this means the pipe flow is pressurized
      result_c[0][c] = PipeD_c[0][c];
    } else if (WettedAngle_c[0][c] < (2.0 * M_PI) && WettedAngle_c[0][c] >= 0.0) {
      // this means the pipe flow is ventilated
      result_c[0][c] = PipeD_c[0][c] * 0.5 * (1.0 - cos(WettedAngle_c[0][c] * 0.5));
    } else {
      // this means the cell is a SW model cell
      result_c[0][c] = WettedArea_c[0][c];
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
WaterDepthEvaluator::EvaluatePartialDerivative_(const State& S,
                                                const Key& wrt_key,
                                                const Tag& wrt_tag,
                                                const std::vector<CompositeVector*>& results)
{
  const auto& WettedAngle_c = *S.Get<CompositeVector>(wetted_angle_key_).ViewComponent("cell");
  const auto& PipeD_c = *S.Get<CompositeVector>(diameter_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == wetted_angle_key_) {
    for (int c = 0; c != ncells; ++c) {
      if (WettedAngle_c[0][c] >= (2.0 * M_PI) || WettedAngle_c[0][c] < 0.0) {
        result_c[0][c] = 0.0;
      } else {
        result_c[0][c] = PipeD_c[0][c] * 0.25 * sin(WettedAngle_c[0][c] * 0.5);
      }
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

} // namespace ShallowWater
} // namespace Amanzi
