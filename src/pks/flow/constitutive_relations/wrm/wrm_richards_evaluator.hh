/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM evaluator calls saturation and rel perm using a capillary pressure p_atm - pc.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_RICHARDS_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_RICHARDS_EVALUATOR_

#include "wrm_evaluator.hh"
#include "wrm.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMRichardsEvaluator : public WRMEvaluator {
 public:

  explicit
  WRMRichardsEvaluator(Teuchos::ParameterList& plist);

  WRMRichardsEvaluator(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<WRMRegionPairList>& wrms);

  WRMRichardsEvaluator(const WRMRichardsEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Key pres_key_;
  bool calc_other_sat_;

};

} // namespae
} // namespae
} // namespae

#endif
