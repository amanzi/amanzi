/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  This WRM model evaluates the saturation of ice, water, and gas.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "wrm_permafrost_model.hh"
#include "secondary_variables_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMPermafrostEvaluator : public SecondaryVariablesFieldEvaluator {
 public:

  explicit
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<WRMPartition>& wrms);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<WRMPermafrostModelPartition>& models);
  WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }
  Teuchos::RCP<WRMPermafrostModelPartition> get_WRMPermafrostModels() { return permafrost_models_; }

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

  Teuchos::RCP<WRMPermafrostModelPartition> permafrost_models_;
  Teuchos::RCP<WRMPartition> wrms_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,WRMPermafrostEvaluator> factory_;

};

} // namespace
} // namespace

#endif
