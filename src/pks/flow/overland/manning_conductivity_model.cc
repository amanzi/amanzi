/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow according to a Manning approach.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "independent_field_model.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ManningConductivityModel::ManningConductivityModel(Teuchos::ParameterList& cond_plist) :
    cond_plist_(cond_plist) {
  slope_regularization_ = cond_plist_.get<double>("slope regularization epsilon", 1.e-8);
  manning_exp_ = cond_plist_.get<double>("Manning exponent", 0.6666666666666667);
  pres_key_ = "overland_pressure";
  slope_key_ = "slope_magnitude";

  my_key_ = "overland_conductivity";
  dependencies_.insert(pres_key_);
  dependencies_.insert(slope_key_);
}


ManningConductivityModel::ManningConductivityModel(const ManningConductivityModel& other) :
    SecondaryVariableFieldModel(other),
    cond_plist_(other.cond_plist_),
    slope_regularization_(other.slope_regularization_),
    pres_key_(other.pres_key_),
    slope_key_(other.slope_key_),
    manning_exp_(other.manning_exp_) {}


Teuchos::RCP<FieldModel>
ManningConductivityModel::Clone() const {
  return Teuchos::rcp(new ManningConductivityModel(*this));
}


// Required methods from SecondaryVariableFieldModel
void ManningConductivityModel::EvaluateField_(const Teuchos::Ptr<State>& S,
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


void ManningConductivityModel::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  ASSERT(0);
  // not implemented, likely not needed.
}


void ManningConductivityModel::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // special EnsureCompatibility to add in a model for Manning's Coef.

  // add the model.
  S->RequireField("manning_coefficient");
  dependencies_.insert("manning_coefficient");
  Teuchos::ParameterList mann_plist = cond_plist_.sublist("Manning coefficient");
  Teuchos::RCP<FieldModel> mann_fm = Teuchos::rcp(new IndependentFieldModel(mann_plist));
  S->SetFieldModel("manning_coefficient", mann_fm);

  // Call the base class's method.
  SecondaryVariableFieldModel::EnsureCompatibility(S);
};


} //namespace
} //namespace
} //namespace

