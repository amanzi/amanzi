/*
  The micropore-macropore flux evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    evalName = micropore_macropore_flux
    modelMethodDeclaration =   double MicroporeMacroporeFlux(double pm, double pM, double krM, double krm, double K) const;
    namespaceCaps = SURFACEBALANCE
    namespace = SurfaceBalance
    evalNameCaps = MICROPORE_MACROPORE_FLUX
    myMethodArgs = pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]
    myKeyMethod = MicroporeMacroporeFlux
    myKeyFirst = micropore
    evalNameString = micropore-macropore flux
    myMethodDeclarationArgs = double pm, double pM, double krM, double krm, double K
    evalClassName = MicroporeMacroporeFlux
    myKey = micropore_macropore_flux  
  
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "micropore_macropore_flux_evaluator.hh"
#include "micropore_macropore_flux_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
MicroporeMacroporeFluxEvaluator::MicroporeMacroporeFluxEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("micropore_macropore_flux parameters");
  model_ = Teuchos::rcp(new MicroporeMacroporeFluxModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
MicroporeMacroporeFluxEvaluator::MicroporeMacroporeFluxEvaluator(const MicroporeMacroporeFluxEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pm_key_(other.pm_key_),
    pM_key_(other.pM_key_),
    krM_key_(other.krM_key_),
    krm_key_(other.krm_key_),
    K_key_(other.K_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
MicroporeMacroporeFluxEvaluator::Clone() const
{
  return Teuchos::rcp(new MicroporeMacroporeFluxEvaluator(*this));
}


// Initialize by setting up dependencies
void
MicroporeMacroporeFluxEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);

  std::string my_key_first("micropore");
  if (domain_name == my_key_first) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // - pull Keys from plist
  // dependency: micropore_pressure
  pm_key_ = plist_.get<std::string>("micropore pressure key",
          domain_name+std::string("micropore_pressure"));
  dependencies_.insert(pm_key_);

  // dependency: pressure
  pM_key_ = plist_.get<std::string>("pressure key",
          domain_name+std::string("pressure"));
  dependencies_.insert(pM_key_);

  // dependency: relative_permeability
  krM_key_ = plist_.get<std::string>("relative permeability key",
          domain_name+std::string("relative_permeability"));
  dependencies_.insert(krM_key_);

  // dependency: micropore_relative_permeability
  krm_key_ = plist_.get<std::string>("micropore relative permeability key",
          domain_name+std::string("micropore_relative_permeability"));
  dependencies_.insert(krm_key_);

  // dependency: micropore_absolute_permeability
  K_key_ = plist_.get<std::string>("micropore absolute permeability key",
          domain_name+std::string("micropore_absolute_permeability"));
  dependencies_.insert(K_key_);
}


void
MicroporeMacroporeFluxEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> pm = S->GetFieldData(pm_key_);
Teuchos::RCP<const CompositeVector> pM = S->GetFieldData(pM_key_);
Teuchos::RCP<const CompositeVector> krM = S->GetFieldData(krM_key_);
Teuchos::RCP<const CompositeVector> krm = S->GetFieldData(krm_key_);
Teuchos::RCP<const CompositeVector> K = S->GetFieldData(K_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
    const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
    const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
    const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
    const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->MicroporeMacroporeFlux(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
    }
  }
}


void
MicroporeMacroporeFluxEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> pm = S->GetFieldData(pm_key_);
Teuchos::RCP<const CompositeVector> pM = S->GetFieldData(pM_key_);
Teuchos::RCP<const CompositeVector> krM = S->GetFieldData(krM_key_);
Teuchos::RCP<const CompositeVector> krm = S->GetFieldData(krm_key_);
Teuchos::RCP<const CompositeVector> K = S->GetFieldData(K_key_);

  if (wrt_key == pm_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMicroporeMacroporeFluxDMicroporePressure(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == pM_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMicroporeMacroporeFluxDPressure(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krM_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMicroporeMacroporeFluxDRelativePermeability(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krm_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMicroporeMacroporeFluxDMicroporeRelativePermeability(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == K_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMicroporeMacroporeFluxDMicroporeAbsolutePermeability(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
