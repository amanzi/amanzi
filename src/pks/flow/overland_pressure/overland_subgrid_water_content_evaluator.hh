/*
  The overland subgrid water content evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
  Subgrid water content.
    
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOW_OVERLAND_SUBGRID_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_OVERLAND_SUBGRID_WATER_CONTENT_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class OverlandSubgridWaterContentModel;

class OverlandSubgridWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  OverlandSubgridWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandSubgridWaterContentEvaluator(const OverlandSubgridWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<OverlandSubgridWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key pres_key_, cv_key_;
  double M_;
  Key delta_max_key_, delta_ex_key_;
  bool bar_;
  double rollover_;
  Teuchos::RCP<OverlandSubgridWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandSubgridWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
