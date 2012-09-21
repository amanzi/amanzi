/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for enthalpy.
----------------------------------------------------------------------------- */


#include "enthalpy_evaluator.hh"

namespace Amanzi {
namespace Energy {

EnthalpyEvaluator::EnthalpyEvaluator(Teuchos::ParameterList& enth_plist) {
  my_key_ = enth_plist.get<std::string>("enthalpy key", "enthalpy_liquid");

  dependencies_.insert(std::string("molar_density_liquid"));
  dependencies_.insert(std::string("internal_energy_liquid"));
  dependencies_.insert(std::string("pressure"));
};

EnthalpyEvaluator::EnthalpyEvaluator(const EnthalpyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other) {};

Teuchos::RCP<FieldEvaluator>
EnthalpyEvaluator::Clone() const {
  return Teuchos::rcp(new EnthalpyEvaluator(*this));
};


void EnthalpyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> u_l = S->GetFieldData("internal_energy_liquid");

  for (int c=0; c!=result->size("cell"); ++c) {
    (*result)("cell",c) = (*u_l)("cell",c) + (*pres)("cell",c)/(*n_l)("cell",c);
  }
};


void EnthalpyEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0);
  // not implemented
};


} //namespace
} //namespace
