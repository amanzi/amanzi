/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow according to a Manning approach.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "independent_variable_field_evaluator.hh"
#include "manning_conductivity_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ManningConductivityEvaluator::ManningConductivityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  slope_regularization_ = plist_.get<double>("slope regularization epsilon", 1.e-8);
  manning_exp_ = plist_.get<double>("Manning exponent", 0.6666666666666667);
  manning_key_ = "manning_coefficient";
  pres_key_ = "overland_pressure";
  slope_key_ = "slope_magnitude";

  my_key_ = "overland_conductivity";
  setLinePrefix(my_key_+std::string(" evaluator"));

  dependencies_.insert(pres_key_);
  dependencies_.insert(slope_key_);
}


ManningConductivityEvaluator::ManningConductivityEvaluator(const ManningConductivityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    slope_regularization_(other.slope_regularization_),
    pres_key_(other.pres_key_),
    slope_key_(other.slope_key_),
    manning_exp_(other.manning_exp_),
    manning_key_(other.manning_key_) {}


Teuchos::RCP<FieldEvaluator>
ManningConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new ManningConductivityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void ManningConductivityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> mann = S->GetFieldData(manning_key_);

  double exponent = manning_exp_ + 1.0;

  for (int c=0; c!=result->size("cell"); ++c) {
    double scaling = (*mann)("cell",c) *
        std::sqrt(std::max((*slope)("cell",c), slope_regularization_));

    if ((*pres)("cell",c) > 0.0) {
      (*result)("cell",c) = std::pow(std::abs((*pres)("cell",c)), exponent) / scaling;
    } else {
      (*result)("cell",c) = 0.0;
    }
  }
}


void ManningConductivityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  ASSERT(0);
  // not implemented, likely not needed.
}


void ManningConductivityEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // special EnsureCompatibility to add in a evaluator for Manning's Coef.

  // add the evaluator.
  S->RequireField(manning_key_);
  dependencies_.insert(manning_key_);
  Teuchos::ParameterList mann_plist = plist_.sublist("Manning coefficient");
  mann_plist.set("evaluator name", manning_key_);
  Teuchos::RCP<FieldEvaluator> mann_fm = Teuchos::rcp(new IndependentVariableFieldEvaluator(mann_plist));
  S->SetFieldEvaluator(manning_key_, mann_fm);

  // Call the base class's method.
  SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
};


} //namespace
} //namespace
} //namespace

