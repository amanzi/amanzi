/*
  Litter drainage rate.

  A simple model based on relaxation from current water content to a saturated water content.
  
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_LITTER_DRAINAGE_EVALUATOR_HH_
#define AMANZI_RELATIONS_LITTER_DRAINAGE_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class LitterDrainageEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  LitterDrainageEvaluator(Teuchos::ParameterList& plist);

  LitterDrainageEvaluator(const LitterDrainageEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key litter_thickness_key_;
  Key litter_wc_key_;
  Key source_key_;
  double source_coef_;
  Key pd_key_;
  double tau_;
  double wc_sat_;
  double n_liq_;
  bool rewetting_;

 private:
  static Amanzi::Utils::RegisteredFactory<FieldEvaluator,LitterDrainageEvaluator> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
