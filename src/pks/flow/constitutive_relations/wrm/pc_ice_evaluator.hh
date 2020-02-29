/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  PCIceEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_PC_ICE_EVALUATOR_HH_
#define AMANZI_RELATIONS_PC_ICE_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class PCIceWater;

class PCIceEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  PCIceEvaluator(Teuchos::ParameterList& plist);
  PCIceEvaluator(const PCIceEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<PCIceWater> get_PCIceWater() { return model_; }

 protected:
  // the actual model
  Teuchos::RCP<PCIceWater> model_;

  // Keys for fields
  // dependencies
  Key temp_key_;
  Key dens_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PCIceEvaluator> factory_;

};

} // namespace
} // namespace

#endif
