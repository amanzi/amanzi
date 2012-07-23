/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Basic interface of a ViscosityModel.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "viscosity_model.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ViscosityModel::ViscosityModel(Teuchos::ParameterList& visc_plist,
        const Teuchos::Ptr<State>& S) :
    visc_plist_(visc_plist) {

  // Process the list for my dependencies.
  // -- temperature
  //      get the name of the temperature field and put it into my list of dependencies
  temp_key_ = visc_plist_.get<string>("temperature key");
  dependencies_.insert(temp_key_);
  //      ensure there is a temperature field in state
  S->RequireField(temp_key_);

  // Process the list for my provided field.
  my_key_ = visc_plist_.get<string>("viscosity key");

  // Ensure there are fields in state, and that this can own them.
  Teuchos::RCP<CompositeVectorFactory> fac_dens = S->RequireField(my_key_, my_key_);
};

ViscosityModel::ViscosityModel(const ViscosityModel& other) :
    SecondaryVariableFieldModel(other),
    visc_plist_(other.visc_plist_),
    temp_key_(other.temp_key_) {}


void ViscosityModel::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector> & result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the model to calculate density.
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int id=0; id!=result->size(*comp); ++id) {
      (*result)(*comp, id) = Viscosity((*temp)(*comp, id));
    }
  }
}


void ViscosityModel::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  // Loop over components, evaluating partial density wrt temp.
  if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      for (int id=0; id!=result->size(*comp); ++id) {
        (*result)(*comp, id) = DViscosityDT((*temp)(*comp, id));
      }
    }
  } else {
    ASSERT(0);
  }
}


} // namespace
} // namespace
} // namespace
