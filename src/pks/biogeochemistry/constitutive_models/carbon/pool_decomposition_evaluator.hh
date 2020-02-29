/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates carbon pool turnover.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_BGCRELATIONS_POOL_DECOMP_HH_
#define AMANZI_BGCRELATIONS_POOL_DECOMP_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

class PoolDecompositionEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  PoolDecompositionEvaluator(Teuchos::ParameterList& plist);
  PoolDecompositionEvaluator(const PoolDecompositionEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key carbon_key_;
  Key decay_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PoolDecompositionEvaluator> fac_;



};

} // namespace
} // namespace
} // namespace

#endif
