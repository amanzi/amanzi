/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow according to a Manning approach.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_

#include "secondary_variable_field_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class ManningConductivityModel : public SecondaryVariableFieldModel {

 public:
  ManningConductivityModel(Teuchos::ParameterList& cond_plist);
  ManningConductivityModel(const ManningConductivityModel& other);
  Teuchos::RCP<FieldModel> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

private:
  Teuchos::ParameterList cond_plist_;
  Key pres_key_;
  Key slope_key_;
  Key manning_key_;
  double slope_regularization_;
  double manning_exp_;
};

} //namespace
} //namespace
} //namespace

#endif

