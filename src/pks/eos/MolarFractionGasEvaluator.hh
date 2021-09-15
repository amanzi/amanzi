/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Determining the molar fraction of a gas component within a gas mixture.
*/

#ifndef AMANZI_EOS_MOLAR_FRACTION_GAS_EVALUATOR_HH_
#define AMANZI_EOS_MOLAR_FRACTION_GAS_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

#include "EOS_SaturatedVaporPressure.hh"
#include "MolarFractionGasEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

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

  Teuchos::RCP<EOS_SaturatedVaporPressure> get_svp_model() { return svp_model_; }

 protected:
  Key temp_key_;

  Teuchos::RCP<EOS_SaturatedVaporPressure> svp_model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator, MolarFractionGasEvaluator> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
