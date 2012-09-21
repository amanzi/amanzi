/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include "dbc.hh"
#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

Utils::RegisteredFactory<WRM,WRMVanGenuchten> WRMVanGenuchten::factory_("van Genuchten");

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMVanGenuchten::WRMVanGenuchten(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

/* ******************************************************************
 * Relative permeability formula: input is capillary pressure pc.
 ****************************************************************** */
double WRMVanGenuchten::k_relative(double pc) {
  if (pc > pc_transition_) {
    double se(pow(1.0 + pow(alpha_*pc,n_),-m_));
    return sqrt(se) * pow( 1.0 - pow( 1.0 - pow(se,1.0/m_),m_), 2);
  } else if (pc <= 0.0) {
    return 1.0;
  } else {
    double se_pct  (pow(1.0 + pow(alpha_*pc_transition_,n_),-m_));
    double f_pct   (sqrt(se_pct) * pow( 1.0 - pow( 1.0 - pow(se_pct,1.0/m_),m_), 2));

    double fab     ((f_pct - 1.0)/pc_transition_);

    double se_pct1 (pow(1.0 + pow(alpha_*(pc_transition_+1.0),n_),-m_));
    double f_pct1  (sqrt(se_pct1) * pow( 1.0 - pow( 1.0 - pow(se_pct1,1.0/m_),m_), 2));


    return 1.0 + pc*pc*fab/pc_transition_ 
      + pc*pc*(pc-pc_transition_)*(f_pct1-f_pct-2.0*fab)/(pc_transition_*pc_transition_);
  }
}


/* ******************************************************************
 * Saturation formula (3.5)-(3.6).
 ****************************************************************** */
double WRMVanGenuchten::saturation(double pc) {
  if (pc > 0.0) {
    return pow(1.0 + pow(alpha_*pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else { 
    return 1.0;
  } 
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double WRMVanGenuchten::d_saturation(double pc) {
  if (pc > 0.0) {
    return -m_*n_ * pow(1.0 + pow(alpha_*pc, n_), -m_-1.0) * pow(alpha_*pc, n_-1) * alpha_ * (1.0 - sr_);
  } else {
    return 0.0;
  }
}

/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
double WRMVanGenuchten::capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  return (pow(pow(se, -1.0/m_) - 1.0, 1/n_)) / alpha_;
}


// set the width of the smoothing interval for the relative permeability 
// function, setting the width to 0.0 will result in the use of the 
// nonsmooth k_rel function

void  WRMVanGenuchten::set_smoothing_interval_width(double pc_transition) {
  ASSERT(pc_transition_ >= 0.0);
  pc_transition_ = pc_transition;
}


void WRMVanGenuchten::InitializeFromPlist_() {
  m_ = plist_.get<double>("van Genuchten m");
  alpha_ = plist_.get<double>("van Genuchten alpha");
  sr_ = plist_.get<double>("van Genuchten residual saturation", 0.0);
  pc_transition_ = plist_.get<double>("van Genuchten smoothing interval width", 0.0);

  n_ = 1.0 / (1.0 - m_);
};

}  // namespace
}  // namespace
}  // namespace
