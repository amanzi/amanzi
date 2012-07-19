/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Basics of an EOS class... likely isn't specific enough for any given EOS,
  and will need other base models.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

EOS::EOS(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S) :
    eos_plist_(eos_plist) {

  // Process the list for my dependencies.
  // -- temperature
  //      get the name of the temperature field and put it into my list of dependencies
  temp_key_ = eos_plist_.get<string>("temperature key");
  dependencies_.insert(temp_key_);
  //      ensure there is a temperature field in state
  S->RequireField(temp_key_);

  // -- same for pressure
  pres_key_ = eos_plist_.get<string>("pressure key");
  dependencies_.insert(pres_key_);
  S->RequireField(pres_key_);

  // Process the list for my provided field.
  my_key_ = eos_plist_.get<string>("density key");

  // Ensure there are fields in state, and that this can own them.
  Teuchos::RCP<CompositeVectorFactory> fac_dens = S->RequireField(my_key_, my_key_);
};


EOS::EOS(const EOS& other) :
    SecondaryVariableFieldModel(other),
    eos_plist_(other.eos_plist_),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_) {}


void EOS::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector> & result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the model to calculate density.
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int id=0; id!=result->size(*comp); ++id) {
      (*result)(*comp, id) = Density((*temp)(*comp, id), (*pres)(*comp, id));
    }
  }
}


void EOS::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  // Loop over components, evaluating partial density wrt pressure or temp.
  if (wrt_key == pres_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      for (int id=0; id!=result->size(*comp); ++id) {
        (*result)(*comp, id) = DDensityDp((*temp)(*comp, id), (*pres)(*comp, id));
      }
    }
  } else if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      for (int id=0; id!=result->size(*comp); ++id) {
        (*result)(*comp, id) = DDensityDT((*temp)(*comp, id), (*pres)(*comp, id));
      }
    }
  } else {
    ASSERT(0);
  }
}


} // namespace
} // namespace
} // namespace
