/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include "soil_thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilThermalConductivityEvaluator::SoilThermalConductivityEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("soil thermal conductivity key",
            "surface-thermal_conductivity");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  AMANZI_ASSERT(plist_.isSublist("soil thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("soil thermal conductivity parameters");

}


SoilThermalConductivityEvaluator::SoilThermalConductivityEvaluator(
      const SoilThermalConductivityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    temperature_key_(other.temperature_key_){}


Teuchos::RCP<FieldEvaluator>
SoilThermalConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new SoilThermalConductivityEvaluator(*this));
}

void SoilThermalConductivityEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  // get mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      // much more efficient to pull out vectors first
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);

      for (int i=0; i!=ncomp; ++i) {

          result_v[0][i] = 2.2; // change this!!!!

      } // i
    }

}


void SoilThermalConductivityEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"SOIL THERMAL CONDUCITIVITY: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
