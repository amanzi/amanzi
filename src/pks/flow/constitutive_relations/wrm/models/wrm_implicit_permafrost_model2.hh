/*
Author: Ethan Coon

Painter's permafrost model.

 */

#ifndef AMANZI_FLOWRELATIONS_WRM_IMPLICIT_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_IMPLICIT_PERMAFROST_MODEL_

#include "wrm_permafrost_model.hh"
#include "wrm_permafrost_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRM;

class WRMImplicitPermafrostModel2 : public WRMPermafrostModel {

 public:
  WRMImplicitPermafrostModel2(Teuchos::ParameterList& plist) :
      WRMPermafrostModel(plist),
      eps_(1.e-12),
      max_it_(100) {}

  // required methods from the base class
  // sats[0] = s_g, sats[1] = s_l, sats[2] = s_i
  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]);
  void saturations(double pc_liq, double pc_ice, double guess, double (&sats)[3]);

  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3]);
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3]);

  // overload version with provided function evaluation
  void dsaturations_dpc_liq(double s_i, double pc_liq, double pc_ice,
          double (&dsats)[3]);
  void  dsaturations_dpc_ice(double s_i, double pc_liq, double pc_ice,
          double (&dsats)[3]);

 private:
  bool saturations_if_saturated_(double pc_liq, double pc_ice, double (&sats)[3]);
  bool saturations_if_above_freezing_(double pc_liq, double pc_ice, double (&sats)[3]);

  double eps_;
  int max_it_;

 private:
  // Functor for saturations()

 private:
  // this Functor gets used within a root-finding algorithm
  class SatIceFunctor_ {
   public:
    SatIceFunctor_(double pc_liq, double pc_ice,
                   const Teuchos::RCP<WRM>& wrm) :
        pc_liq_(pc_liq), pc_ice_(pc_ice), wrm_(wrm) {}

    double operator()(double s_i) {
      double tmp = (1.0 - s_i) * wrm_->saturation(pc_liq_);
      return tmp - wrm_->saturation( pc_ice_ + wrm_->capillaryPressure( tmp + s_i));
    }

   private:
    double pc_liq_;
    double pc_ice_;
    Teuchos::RCP<WRM> wrm_;

  };

  struct Tol_ {
    Tol_(double eps) : eps_(eps) {}
    bool operator()(const double& a, const double& b) const {
      return std::abs(a - b) <= eps_;
    }
    double eps_;
  };

 private:
  static Utils::RegisteredFactory<WRMPermafrostModel,WRMImplicitPermafrostModel2> factory_;

};


} //namespace
} //namespace
} //namespace

#endif
