/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Model simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_WRM_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_MODEL_

#include "secondary_variables_field_model.hh"
#include "wrm.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMModel : public SecondaryVariablesFieldModel {

 public:
  // constructor format for all derived classes
  WRMModel(Teuchos::ParameterList& wrm_plist);
  WRMModel(const WRMModel& other);

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) = 0;

 protected:

  Teuchos::ParameterList wrm_plist_;
  Teuchos::RCP<WRM> wrm_;
};

} //namespace
} //namespace
} //namespace

#endif
