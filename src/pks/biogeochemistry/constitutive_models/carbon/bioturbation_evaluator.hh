/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates bioturbation of carbon -- simple diffusion model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_BGCRELATIONS_BIOTURBATION_HH_
#define AMANZI_BGCRELATIONS_BIOTURBATION_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

class BioturbationEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  BioturbationEvaluator(Teuchos::ParameterList& plist);
  BioturbationEvaluator(const BioturbationEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

protected:
  Key carbon_key_;
  Key diffusivity_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,BioturbationEvaluator> fac_;



};

} // namespace
} // namespace
} // namespace

#endif
