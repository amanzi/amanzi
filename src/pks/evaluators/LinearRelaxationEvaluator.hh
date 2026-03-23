/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  The linear relaxion evaluator a(t) (u - b(t)), where u is 
  dependent variable.
*/

#ifndef AMANZI_EVALUATORS_LINEAR_RELAXATION_HH_
#define AMANZI_EVALUATORS_LINEAR_RELAXATION_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "Key.hh"
#include "State.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace Evaluators {

class LinearRelaxationEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit LinearRelaxationEvaluator(Teuchos::ParameterList& plist, Teuchos::RCP<State>& S);
  LinearRelaxationEvaluator(const LinearRelaxationEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  std::vector<Teuchos::RCP<PK_DomainFunction>> srcs_;
  Key variable_key_; 
};

}  // namespace Evaluators
}  // namespace Amanzi

#endif

