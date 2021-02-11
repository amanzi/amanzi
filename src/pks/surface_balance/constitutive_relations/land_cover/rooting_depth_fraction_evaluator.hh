/*
  The rooting depth fraction evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.


    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_ROOTING_DEPTH_FRACTION_EVALUATOR_HH_
#define AMANZI_FLOW_ROOTING_DEPTH_FRACTION_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace LandCover {
namespace Relations {

class RootingDepthFractionModel;

class RootingDepthFractionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  RootingDepthFractionEvaluator(Teuchos::ParameterList& plist);
  RootingDepthFractionEvaluator(const RootingDepthFractionEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  std::vector<Teuchos::RCP<RootingDepthFractionModel>> get_models() { return models_; }

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  void InitializeFromPlist_();

  Key z_key_;
  Key cv_key_;
  Key surf_cv_key_;
  int npfts_;

  std::vector<Teuchos::RCP<RootingDepthFractionModel>> models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RootingDepthFractionEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
