/*
  The surface ice energy evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Energy evaulator for ice+liquid surface water.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_SURFACE_ICE_ENERGY_EVALUATOR_HH_
#define AMANZI_ENERGY_SURFACE_ICE_ENERGY_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class SurfaceIceEnergyModel;

class SurfaceIceEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SurfaceIceEnergyEvaluator(Teuchos::ParameterList& plist);
  SurfaceIceEnergyEvaluator(const SurfaceIceEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<SurfaceIceEnergyModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key h_key_;
  Key eta_key_;
  Key nl_key_;
  Key ul_key_;
  Key ni_key_;
  Key ui_key_;
  Key cv_key_;

  Teuchos::RCP<SurfaceIceEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SurfaceIceEnergyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif