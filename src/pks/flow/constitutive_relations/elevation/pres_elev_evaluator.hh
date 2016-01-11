/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_PRES_ELEV_EVALUATOR_
#define AMANZI_FLOWRELATIONS_PRES_ELEV_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class PresElevEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PresElevEvaluator(Teuchos::ParameterList& plist);

  PresElevEvaluator(const PresElevEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;
  
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  Key pres_key_;
  Key elev_key_;
};

} //namespace
} //namespace
} //namespace

#endif
