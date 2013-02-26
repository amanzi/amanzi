/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity after a cell volume change due to mesh deformation.

  Authors: Markus Berndt (berndt@lanl.gov)
*/

#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {
namespace DeformRelations {


Utils::RegisteredFactory<FieldEvaluator,PorosityEvaluator> PorosityEvaluator::factory_("deformation porosity");

PorosityEvaluator::PorosityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  my_key_ = "porosity";
  setLinePrefix(my_key_+std::string(" evaluator"));
  
  // add dependency to cell volume
  dependencies_.insert("cell_volume");
  dependencies_.insert("deformation");
}


PorosityEvaluator::PorosityEvaluator(const PorosityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other)
{ }

Teuchos::RCP<FieldEvaluator>
PorosityEvaluator::Clone() const {
  return Teuchos::rcp(new PorosityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void PorosityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  CompositeVector& rho = *S->GetFieldData("porosity",my_key_);
  const CompositeVector& deformation = *S->GetFieldData("deformation");

  // (new rho) = deformation + (old rho) - 1.0
  rho.Update(1.0, deformation, 1.0);
  rho.Shift(-1.0);
}



void PorosityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  ASSERT(0);
  // not implemented, likely not needed.
}

  
void PorosityEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // special EnsureCompatibility to add in a evaluator for Porosity
  
  // Call the base class's method since we do not need anything special here
  SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
};


} //namespace
} //namespace
} //namespace

