/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_VISC_EVALUATOR_HH_
#define AMANZI_RELATIONS_VISC_EVALUATOR_HH_

#include "viscosity_relation.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class ViscosityEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  ViscosityEvaluator(Teuchos::ParameterList& plist);

  ViscosityEvaluator(const ViscosityEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // the actual model
  Teuchos::RCP<ViscosityRelation> visc_;

  // Keys for fields
  // dependencies
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ViscosityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
