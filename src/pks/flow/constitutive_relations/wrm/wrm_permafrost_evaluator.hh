/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM model evaluates the saturation of ice, water, and gas.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_

#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMPermafrostEvaluator : public SecondaryVariablesFieldEvaluator {
 public:

  explicit
  WRMPermafrostEvaluator(Teuchos::ParameterList& wrm_plist);

  WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 private:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 private:
  Teuchos::ParameterList wrm_plist_;
  Key one_on_A_key_;
  Key one_on_B_key_;
  Key s_l_key_;

};

} // namespace
} // namespace
} // namespace

#endif
