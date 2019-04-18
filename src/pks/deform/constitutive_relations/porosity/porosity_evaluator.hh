/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow according to a Manning approach.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_DEFORMRELATIONS_POROSITY_EVALUATOR_
#define AMANZI_DEFORMRELATIONS_POROSITY_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"



namespace Amanzi {
namespace Deform {
namespace DeformRelations {

class PorosityEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PorosityEvaluator(Teuchos::ParameterList& cond_plist);
  PorosityEvaluator(const PorosityEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

private:
  static Utils::RegisteredFactory<FieldEvaluator,PorosityEvaluator> factory_;


};

} //namespace
} //namespace
} //namespace

#endif

