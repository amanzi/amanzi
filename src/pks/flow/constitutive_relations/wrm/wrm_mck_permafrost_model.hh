/*
Author: Ethan Coon

McKenzie et al. (2007)'s soil freezing curve

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_MCK_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_MCK_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {

class WRM;

class WRMMCKPermafrostModel : public WRMPermafrostModel {

 public:
  double residualSaturation();
  explicit
  WRMMCKPermafrostModel(Teuchos::ParameterList& plist) :
      WRMPermafrostModel(plist) {
    T0_ = plist.get<double>("freezing point [K]", 273.15);
  //  sr_ = plist.get<double>("residual saturation [-]");
    w_  = plist.get<double>("sfc fitting coefficient", 3.0);
  //  AMANZI_ASSERT(sr_ >= 0.); AMANZI_ASSERT(sr_ < 1.);
  }


  // required methods from the base class
  // sats[0] = sg, sats[1] = sl, sats[2] = si
  virtual bool freezing(double T, double pc_liq, double pc_ice) { 
    return T < T0_;
  }

  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

 protected:
  double T0_;
  double w_;
  double sr_;

 private:
  // factory registration
  static Utils::RegisteredFactory<WRMPermafrostModel,WRMMCKPermafrostModel> factory_;

};


} //namespace
} //namespace

#endif
