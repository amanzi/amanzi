/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {}

MultiplicativeEvaluator::MultiplicativeEvaluator(const MultiplicativeEvaluator& other) :
    SecondaryVariableFieldEvaluator(other) {}

Teuchos::RCP<FieldEvaluator>
MultiplicativeEvaluator::Clone() const
{
  return Teuchos::rcp(new MultiplicativeEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void
MultiplicativeEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
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

void
MultiplicativeEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
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


} // namespace
} // namespace

