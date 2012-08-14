/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM model calls saturation and rel perm using a capillary pressure p_atm - pc.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_WRM_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_PERMAFROST_MODEL_

#include "wrm_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMPermafrostModel : public SecondaryVariablesFieldModel {
 public:

  explicit
  WRMPermafrostModel(Teuchos::ParameterList& wrm_plist);

  WRMPermafrostModel(const WRMPermafrostModel& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

 private:
  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 private:
  Teuchos::ParameterList wrm_plist_;
  Key one_on_A_key_;
  Key one_on_B_key_;
  Key s_l_key_;

};

} // namespae
} // namespae
} // namespae

#endif
