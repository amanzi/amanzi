/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM model calls saturation using a thermal term.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_WRM_ICE_WATER_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_ICE_WATER_MODEL_

#include "pc_ice_water.hh"
#include "wrm_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRM; // forward declaration

class WRMIceLiquidModel : public WRMModel {
 public:

  explicit
  WRMIceLiquidModel(Teuchos::ParameterList& wrm_plist);

  WRMIceLiquidModel(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm);

  WRMIceLiquidModel(const WRMIceLiquidModel& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

 protected:
  Key temp_key_;
  Key dens_key_;
  bool calc_other_sat_;

  Teuchos::RCP<PCIceWater> pc_;
};

} // namespae
} // namespae
} // namespae

#endif
