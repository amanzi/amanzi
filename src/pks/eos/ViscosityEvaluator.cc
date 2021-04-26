/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  ViscosityEvaluator is the interface between state/data and the 
  model, a VPM.
*/

#include "ViscosityBaseFactory.hh"
#include "ViscosityEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

ViscosityEvaluator::ViscosityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  // my keys
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("viscosity key", "viscosity_liquid");
  }

  // set up my dependencies
  std::string domain = Keys::getDomain(my_key_);
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(temp_key_);

  pres_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pres_key_);

  // Construct my Viscosity model
  AMANZI_ASSERT(plist_.isSublist("viscosity model parameters"));
  ViscosityBaseFactory visc_fac;
  visc_ = visc_fac.CreateViscosity(plist_.sublist("viscosity model parameters"));
}


ViscosityEvaluator::ViscosityEvaluator(const ViscosityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    visc_(other.visc_),
    temp_key_(other.temp_key_) {};


Teuchos::RCP<FieldEvaluator> ViscosityEvaluator::Clone() const {
  return Teuchos::rcp(new ViscosityEvaluator(*this));
}


void ViscosityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result)
{
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // evaluate p_s / p_atm
  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp));
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      AMANZI_ASSERT(temp_v[0][id] > 200.);
      result_v[0][id] = visc_->Viscosity(temp_v[0][id], pres_v[0][id]);
    }
  }
}


void ViscosityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(wrt_key == temp_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // evaluate d/dT( p_s / p_atm )
  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp));
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = visc_->DViscosityDT(temp_v[0][id], pres_v[0][id]);
    }
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi
