/*
  The interception fraction evaluator is an algebraic evaluator of a given model.
  
  Generated via evaluator_generator.
*/

#include "interception_fraction_evaluator.hh"
#include "interception_fraction_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
InterceptionFractionEvaluator::InterceptionFractionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("interception_fraction parameters");
  model_ = Teuchos::rcp(new InterceptionFractionModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
InterceptionFractionEvaluator::InterceptionFractionEvaluator(const InterceptionFractionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    ai_key_(other.ai_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
InterceptionFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new InterceptionFractionEvaluator(*this));
}


// Initialize by setting up dependencies
void
InterceptionFractionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomainPrefix(my_key_);

  // - pull Keys from plist
  // dependency: surface-area_index
  ai_key_ = plist_.get<std::string>("area index key", domain_name+"area_index");
  dependencies_.insert(ai_key_);
}


void
InterceptionFractionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> ai = S->GetFieldData(ai_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& ai_v = *ai->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->InterceptionFraction(ai_v[0][i]);
    }
  }
}


void
InterceptionFractionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> ai = S->GetFieldData(ai_key_);

  if (wrt_key == ai_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& ai_v = *ai->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DInterceptionFractionDAreaIndex(ai_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
