/*
  The evaporation downregulation via soil resistance evaluator is an algebraic evaluator of a given model.
Downregulates evaporation from a potential.
  
  Generated via evaluator_generator.
*/

#include "evaporation_downregulation_evaluator.hh"
#include "evaporation_downregulation_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporationDownregulationEvaluator::EvaporationDownregulationEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("evaporation_downregulation parameters");
  model_ = Teuchos::rcp(new EvaporationDownregulationModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
EvaporationDownregulationEvaluator::EvaporationDownregulationEvaluator(const EvaporationDownregulationEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    sg_key_(other.sg_key_),
    poro_key_(other.poro_key_),
    pot_evap_key_(other.pot_evap_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
EvaporationDownregulationEvaluator::Clone() const
{
  return Teuchos::rcp(new EvaporationDownregulationEvaluator(*this));
}


// Initialize by setting up dependencies
void
EvaporationDownregulationEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomainPrefix(my_key_);

  // - pull Keys from plist
  // dependency: saturation_gas
  sg_key_ = plist_.get<std::string>("saturation gas key",
          domain_name+"saturation_gas");
  dependencies_.insert(sg_key_);

  // dependency: porosity
  poro_key_ = plist_.get<std::string>("porosity key",
          domain_name+"porosity");
  dependencies_.insert(poro_key_);

  // dependency: potential_evaporation
  pot_evap_key_ = plist_.get<std::string>("potential evaporation key",
          domain_name+"potential_evaporation");
  dependencies_.insert(pot_evap_key_);
}


void
EvaporationDownregulationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> sg = S->GetFieldData(sg_key_);
Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
Teuchos::RCP<const CompositeVector> pot_evap = S->GetFieldData(pot_evap_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp, false);
    const Epetra_MultiVector& pot_evap_v = *pot_evap->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Evaporation(sg_v[0][i], poro_v[0][i], pot_evap_v[0][i]);
    }
  }
}


void
EvaporationDownregulationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> sg = S->GetFieldData(sg_key_);
Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
Teuchos::RCP<const CompositeVector> pot_evap = S->GetFieldData(pot_evap_key_);

  if (wrt_key == sg_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp, false);
      const Epetra_MultiVector& pot_evap_v = *pot_evap->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEvaporationDSaturationGas(sg_v[0][i], poro_v[0][i], pot_evap_v[0][i]);
      }
    }

  } else if (wrt_key == poro_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp, false);
      const Epetra_MultiVector& pot_evap_v = *pot_evap->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEvaporationDPorosity(sg_v[0][i], poro_v[0][i], pot_evap_v[0][i]);
      }
    }

  } else if (wrt_key == pot_evap_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp, false);
      const Epetra_MultiVector& pot_evap_v = *pot_evap->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEvaporationDPotentialEvaporation(sg_v[0][i], poro_v[0][i], pot_evap_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
