/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_ICY_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_ICY_HEIGHT_EVALUATOR_

#include "height_evaluator.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class IcyHeightEvaluator : public HeightEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  IcyHeightEvaluator(Teuchos::ParameterList& plist);
  IcyHeightEvaluator(const IcyHeightEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key dens_ice_key_;
  Key unfrozen_frac_key_;

};

} //namespace
} //namespace
} //namespace

#endif
