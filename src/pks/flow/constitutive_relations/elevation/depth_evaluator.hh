/*
  Evaluates depth of various mesh entities.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_DEPTH_EVALUATOR_HH_
#define AMANZI_FLOW_DEPTH_EVALUATOR_HH_

#include "Factory.hh"
#include "independent_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class DepthEvaluator : public IndependentVariableFieldEvaluator {

 public:
  explicit
  DepthEvaluator(Teuchos::ParameterList& plist);
  DepthEvaluator(const DepthEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

 protected:
  // Required methods from IndependentVariableFieldEvaluator
  virtual void UpdateField_(const Teuchos::Ptr<State>& S) override;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,DepthEvaluator> reg_;

};

} //namespace
} //namespace

#endif
