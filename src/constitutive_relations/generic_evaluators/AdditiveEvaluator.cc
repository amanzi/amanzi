/*
  AdditiveEvaluator is the generic evaluator for adding N other fields.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

AdditiveEvaluator::AdditiveEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::Array<std::string> deps =
      plist_.get<Teuchos::Array<std::string> >("evaluator dependencies");

  for (Teuchos::Array<std::string>::iterator dep=deps.begin();
       dep!=deps.end(); ++dep) {
    Key pname = *dep + std::string(" coefficient");
    coefs_[*dep] = plist.get<double>(pname, 1.0);
  }  
}


AdditiveEvaluator::AdditiveEvaluator(const AdditiveEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    coefs_(other.coefs_) {}

Teuchos::RCP<FieldEvaluator>
AdditiveEvaluator::Clone() const
{
  return Teuchos::rcp(new AdditiveEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
AdditiveEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  ASSERT(dependencies_.size() > 1);
  result->PutScalar(0.);
  
  for (std::map<Key, double>::const_iterator it=coefs_.begin();
       it!=coefs_.end(); ++it) {
    Teuchos::RCP<const CompositeVector> dep = S->GetFieldData(it->first);
    result->Update(it->second, *dep, 1.0);
  }
}

void
AdditiveEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(coefs_[wrt_key]);
}


} // namespace
} // namespace

