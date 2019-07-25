/*
  The transpiration distribution via rooting depth evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Distributes transpiration based upon a rooting depth and a wilting-point water-potential factor.

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_TRANSPIRATION_DISTRIBUTION_EVALUATOR_HH_
#define AMANZI_FLOW_TRANSPIRATION_DISTRIBUTION_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

class Function;

namespace Flow {
namespace Relations {

class TranspirationDistributionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  TranspirationDistributionEvaluator(Teuchos::ParameterList& plist);
  TranspirationDistributionEvaluator(const TranspirationDistributionEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
 protected:
  void InitializeFromPlist_();

  Key f_wp_key_;
  Key f_root_key_;
  Key trans_total_key_;
  Key cv_key_;
  Key surf_cv_key_;
  int npfts_;

  bool limiter_local_;
  Teuchos::RCP<Function> limiter_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,TranspirationDistributionEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
