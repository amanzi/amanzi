/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Amanzi needs to fix its WRM interface... the current lack of factory or
  common interface is stupid.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_WRM_VAN_GENUCHTEN_
#define _FLOWRELATIONS_WRM_VAN_GENUCHTEN_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMVanGenuchten : public WRM {

public:
  explicit WRMVanGenuchten(Teuchos::ParameterList& wrm_plist);

  // required methods from the base class
  double k_relative(double pc);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

  void set_smoothing_interval_width(double width);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList wrm_plist_;

  double m_;  // van Genuchten parameters: m, n, alpha
  double n_;
  double alpha_;
  double sr_;  // van Genuchten residual saturation

  // the following is for smoothing the saturation and k_relative curves
  double pc_transition_; // we smooth curves in the interval [0, pc_transition]

  static Utils::RegisteredFactory<WRM,WRMVanGenuchten> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
