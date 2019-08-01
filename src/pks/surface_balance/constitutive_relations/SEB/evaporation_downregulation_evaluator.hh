/*
  The evaporation downregulation via soil resistance evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Downregulates evaporation from a potential.

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SEB_EVAPORATION_DOWNREGULATION_EVALUATOR_HH_
#define AMANZI_SEB_EVAPORATION_DOWNREGULATION_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporationDownregulationModel;

class EvaporationDownregulationEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  EvaporationDownregulationEvaluator(Teuchos::ParameterList& plist);
  EvaporationDownregulationEvaluator(const EvaporationDownregulationEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<EvaporationDownregulationModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key sg_key_;
  Key poro_key_;
  Key pot_evap_key_;

  Teuchos::RCP<EvaporationDownregulationModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EvaporationDownregulationEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif