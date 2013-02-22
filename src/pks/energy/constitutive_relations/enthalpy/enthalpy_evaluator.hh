/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for enthalpy.
----------------------------------------------------------------------------- */


#ifndef AMANZI_ENTHALPY_EVALUATOR_HH_
#define AMANZI_ENTHALPY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {

class EnthalpyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  EnthalpyEvaluator(Teuchos::ParameterList& plist);
  EnthalpyEvaluator(const EnthalpyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key pres_key_;
  Key dens_key_;
  Key ie_key_;
  bool include_work_;
};

} // namespace
} // namespace

#endif
