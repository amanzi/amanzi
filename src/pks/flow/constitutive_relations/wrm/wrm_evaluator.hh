/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_

#include "wrm.hh"
#include "secondary_variables_field_evaluator.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  WRMEvaluator(Teuchos::ParameterList& plist);
  WRMEvaluator(Teuchos::ParameterList& plist,
               const Teuchos::RCP<WRMRegionPairList>& wrms);
  WRMEvaluator(const WRMEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;


 protected:
  void InitializeFromPlist_();
  Teuchos::RCP<WRMRegionPairList> get_WRMs() { return wrms_; }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Teuchos::RCP<WRMRegionPairList> wrms_;
  bool calc_other_sat_;
  Key cap_pres_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> factory_;

};

} //namespace
} //namespace
} //namespace

#endif
