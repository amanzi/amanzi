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

ViscosityModel::ViscosityModel(Teuchos::ParameterList& visc_plist) :
    visc_plist_(visc_plist) {

  // Process the list for my provided field.
  // Process the list for my provided field.
  if (visc_plist_.isParameter("viscosity key")) {
    my_key_ = visc_plist_.get<string>("viscosity key");
  } else {
    std::string name = visc_plist_.name();
    std::size_t start = name.find_last_of(">");
    my_key_ = name.substr(start+1);
  }

  // Set up my dependencies.
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);
  if (domain_name == std::string("viscosity")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // -- temperature
  temp_key_ = domain_name+std::string("temperature");
  dependencies_.insert(temp_key_);
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
  // evaluating the model to calculate viscosity.
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

  // Loop over components, evaluating partial viscosity wrt temp.
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
