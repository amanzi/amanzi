/*
Author: Ethan Coon

Painter's permafrost model.

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_OLD_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_OLD_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRM;

class WRMOldPermafrostModel : public WRMPermafrostModel {

 public:
  WRMOldPermafrostModel(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<WRM>& wrm) :
      plist_(plist),
      wrm_(wrm) {}

  // required methods from the base class
  // sats[0] = s_g, sats[1] = s_l, sats[2] = s_i
  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 private:

  Teuchos::ParameterList plist_;
  Teuchos::RCP<WRM> wrm_;
};


} //namespace
} //namespace
} //namespace

#endif
