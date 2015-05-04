/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_MULTIPLICATIVE_EVALUATOR_
#define AMANZI_RELATIONS_MULTIPLICATIVE_EVALUATOR_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class MultiplicativeEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
      SecondaryVariableFieldEvaluator(plist) {}

  MultiplicativeEvaluator(const MultiplicativeEvaluator& other) :
      SecondaryVariableFieldEvaluator(other) {}

  Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new MultiplicativeEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  void EvaluateField_(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<CompositeVector>& result) {
    ASSERT(dependencies_.size() > 1);
    KeySet::const_iterator key = dependencies_.begin();
    *result = *S->GetFieldData(*key);
    key++;

    for (; key!=dependencies_.end(); ++key) {
      Teuchos::RCP<const CompositeVector> dep = S->GetFieldData(*key);
      for (CompositeVector::name_iterator lcv=result->begin(); lcv!=result->end(); ++lcv) {
        Epetra_MultiVector& res_c = *result->ViewComponent(*lcv, false);
        const Epetra_MultiVector& dep_c = *dep->ViewComponent(*lcv, false);

        for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
      }      
    }
  }

  void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    ASSERT(dependencies_.size() > 1);

    KeySet::const_iterator key = dependencies_.begin();
    while (*key == wrt_key) key++;
    *result = *S->GetFieldData(*key);
    key++;

    for (; key!=dependencies_.end(); ++key) {
      if (*key != wrt_key) {
        Teuchos::RCP<const CompositeVector> dep = S->GetFieldData(*key);
        for (CompositeVector::name_iterator lcv=result->begin(); lcv!=result->end(); ++lcv) {
          Epetra_MultiVector& res_c = *result->ViewComponent(*lcv, false);
          const Epetra_MultiVector& dep_c = *dep->ViewComponent(*lcv, false);

          for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
        }      
      }
    }
  }
    
 private:
  static Utils::RegisteredFactory<FieldEvaluator,MultiplicativeEvaluator> factory_;
};

} // namespace
} // namespace

#endif
