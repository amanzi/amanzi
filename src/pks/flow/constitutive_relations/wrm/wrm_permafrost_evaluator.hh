/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM model evaluates the saturation of ice, water, and gas.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_

#include "wrm.hh"
#include "wrm_permafrost_model.hh"
#include "secondary_variables_field_evaluator.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMPermafrostEvaluator : public SecondaryVariablesFieldEvaluator {
 public:

  explicit
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<WRMRegionPairList>& wrms);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<WRMPermafrostModelRegionPairList>& models);
  WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

  void InitializeFromPlist_();

 protected:
  Key pc_liq_key_;
  Key pc_ice_key_;
  Key s_l_key_;

  Teuchos::RCP<WRMPermafrostModelRegionPairList> permafrost_models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,WRMPermafrostEvaluator> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
