/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.
*/

#ifndef AMANZI_EOS_EFFECTIVE_PRESSURE_EVALUATOR_HH_
#define AMANZI_EOS_EFFECTIVE_PRESSURE_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class EffectivePressureEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  explicit
  EffectivePressureEvaluator(Teuchos::ParameterList& ep_plist);

  EffectivePressureEvaluator(const EffectivePressureEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Teuchos::ParameterList ep_plist_;

  // Keys for fields dependencies
  Key pres_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EffectivePressureEvaluator> factory_;
};

}  // namespace Relations
}  // namespace Amanzi

#endif
