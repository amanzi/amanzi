/*
  The rooting depth fraction evaluator is an algebraic evaluator of a given model.
Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

  
  Generated via evaluator_generator.
*/

#include "rooting_depth_fraction_evaluator.hh"
#include "rooting_depth_fraction_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("rooting_depth_fraction parameters");
  model_ = Teuchos::rcp(new RootingDepthFractionModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(const RootingDepthFractionEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    z_key_(other.z_key_),    
    cv_key_(other.cv_key_),    
    surface_cv_key_(other.surface_cv_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
RootingDepthFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new RootingDepthFractionEvaluator(*this));
}


// Initialize by setting up dependencies
void
RootingDepthFractionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomainPrefix(my_key_);

  // - pull Keys from plist
  // dependency: depth
  z_key_ = Keys::readKey(plist_, domain_name, "depth", "depth");
  dependencies_.insert(z_key_);

  // Manually added:
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");

  Key surf_domain_name;
  if (domain_name.empty() || domain_name == "domain") {
    surf_domain_name = "surface";
  } else {
    surf_domain_name = Key("surface_")+domain_name;
  }
  surf_domain_name = plist_.get<std::string>("surface domain name", surf_domain_name);
    
  surface_cv_key_ = Keys::readKey(plist_, surf_domain_name, "cell volume", "cell_volume");
}


void
RootingDepthFractionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> z = S->GetFieldData(z_key_);
  Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(cv_key_);
  Teuchos::RCP<const CompositeVector> surf_cv = S->GetFieldData(surface_cv_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    ASSERT(*comp == "cell");
    const Epetra_MultiVector& z_v = *z->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);
    Epetra_MultiVector integral(result_v);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->RootingDepthFraction(z_v[0][i]);
    }
  }
}


void
RootingDepthFractionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> z = S->GetFieldData(z_key_);

  if (wrt_key == z_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& z_v = *z->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DRootingDepthFractionDDepth(z_v[0][i]);
      }
    }

  } else {
    ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
