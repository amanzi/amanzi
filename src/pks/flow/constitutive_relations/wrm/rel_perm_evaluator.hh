/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Rel perm( pc ( sat ) ).

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_
#define AMANZI_FLOWRELATIONS_REL_PERM_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "wrm.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class RelPermEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  //  explicit
  //  RelPermEvaluator(Teuchos::ParameterList& plist);
  RelPermEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<WRMRegionPairList>& wrms);

  RelPermEvaluator(const RelPermEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<WRMRegionPairList> get_WRMs() { return wrms_; }

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMRegionPairList> wrms_;
  Key sat_key_;

  double min_val_;
};

} //namespace
} //namespace
} //namespace

#endif
