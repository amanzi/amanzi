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

#ifndef AMANZI_WATER_DEPTH_EVALUATOR_HH_
#define AMANZI_WATER_DEPTH_EVALUATOR_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ShallowWater {

class WaterDepthEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  WaterDepthEvaluator(Teuchos::ParameterList& plist);

  // required interface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  std::string wetted_angle_key_;
  std::string primary_variable_key_;
  std::string diameter_key_;
};

} // namespace ShallowWater
} // namespace Amanzi

#endif
