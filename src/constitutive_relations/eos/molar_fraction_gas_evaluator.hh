/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONSRELATIONS_MOLAR_FRACTION_GAS_
#define AMANZI_RELATIONSRELATIONS_MOLAR_FRACTION_GAS_

#include "vapor_pressure_relation.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class MolarFractionGasEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  MolarFractionGasEvaluator(Teuchos::ParameterList& plist);

  MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key temp_key_;

  Teuchos::RCP<VaporPressureRelation> sat_vapor_model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,MolarFractionGasEvaluator> factory_;

};

} //namespace
} //namespace

#endif
