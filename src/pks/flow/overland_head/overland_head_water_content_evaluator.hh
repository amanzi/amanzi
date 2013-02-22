/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_WATER_CONTENT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_OVERLAND_HEAD_WATER_CONTENT_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class OverlandHeadWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  OverlandHeadWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandHeadWaterContentEvaluator(const OverlandHeadWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key dens_key_;
  Key height_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandHeadWaterContentEvaluator> factory_;

};

} //namespace
} //namespace
} //namespace

#endif
