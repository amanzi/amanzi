/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Determining the molar fraction of a gas component within a gas mixture.
*/

#ifndef AMANZI_EOS_MOLAR_FRACTION_GAS_HH_
#define AMANZI_EOS_MOLAR_FRACTION_GAS_HH_

#include "vapor_pressure_relation.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class MolarFractionGasEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit MolarFractionGasEvaluator(Teuchos::ParameterList& plist);

  MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<VaporPressureRelation> get_VaporPressureRelation() {
    return sat_vapor_model_; }

 protected:
  Key temp_key_;

  Teuchos::RCP<VaporPressureRelation> sat_vapor_model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,MolarFractionGasEvaluator> factory_;
};

}  // namespace Relations
}  // namespace Amanzi

#endif
