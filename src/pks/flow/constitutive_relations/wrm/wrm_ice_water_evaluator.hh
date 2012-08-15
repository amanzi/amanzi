/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  This WRM model calls the wrm model using a capillary pressure given by first
  calculating the capillary pressure of water on ice using a thermal relation
  and then shifting it to the gas on water capillary pressure using a ratio of
  surface tensions.  See the permafrost notes.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_WRM_ICE_WATER_EVALUATOR_
#define AMANZI_FLOWRELATIONS_WRM_ICE_WATER_EVALUATOR_

#include "pc_ice_water.hh"
#include "wrm.hh"
#include "wrm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMIceWaterEvaluator : public WRMEvaluator {
 public:

  explicit
  WRMIceWaterEvaluator(Teuchos::ParameterList& wrm_plist);

  WRMIceWaterEvaluator(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm);

  WRMIceWaterEvaluator(const WRMIceWaterEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
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
