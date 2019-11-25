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

#ifndef AMANZI_FLOW_ONE_PFT_ROOTING_DEPTH_FRACTION_EVALUATOR_HH_
#define AMANZI_FLOW_ONE_PFT_ROOTING_DEPTH_FRACTION_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class RootingDepthFractionModel;

class OnePFTRootingDepthFractionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  OnePFTRootingDepthFractionEvaluator(Teuchos::ParameterList& plist);
  OnePFTRootingDepthFractionEvaluator(const OnePFTRootingDepthFractionEvaluator& other) = default;
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  std::vector<std::pair<std::string, Teuchos::RCP<RootingDepthFractionModel> > > get_models() { return models_; }

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  void InitializeFromPlist_();

  Key z_key_;
  Key cv_key_;
  Key surf_cv_key_;

  Key surf_domain_;
  Key subsurf_domain_;

  std::vector<std::pair<std::string, Teuchos::RCP<RootingDepthFractionModel> > > models_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OnePFTRootingDepthFractionEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
