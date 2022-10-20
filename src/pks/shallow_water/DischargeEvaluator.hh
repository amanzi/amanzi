/*
  Shallow Water PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DISCHARGE_EVALUATOR_HH_
#define AMANZI_DISCHARGE_EVALUATOR_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ShallowWater {

class DischargeEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  DischargeEvaluator(Teuchos::ParameterList& plist);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  std::string velocity_key_, ponded_depth_key_;
};

}  // namespace ShallowWater
}  // namespace Amanzi

#endif

