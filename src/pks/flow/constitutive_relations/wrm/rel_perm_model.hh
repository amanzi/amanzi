/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Rel perm of sat_l.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_REL_PERM_MODEL_
#define AMANZI_FLOWRELATIONS_REL_PERM_MODEL_

#include "secondary_variable_field_model.hh"
#include "wrm.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class RelPermModel : public SecondaryVariableFieldModel {

 public:
  // constructor format for all derived classes
  explicit
  RelPermModel(Teuchos::ParameterList& wrm_plist);
  RelPermModel(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm);
  RelPermModel(const RelPermModel& other);
  virtual Teuchos::RCP<FieldModel> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<WRM> get_WRM() { return wrm_; }

 protected:
  void InitializeFromPlist_();

  Teuchos::ParameterList wrm_plist_;
  Teuchos::RCP<WRM> wrm_;
  Key sat_key_;
};

} //namespace
} //namespace
} //namespace

#endif
