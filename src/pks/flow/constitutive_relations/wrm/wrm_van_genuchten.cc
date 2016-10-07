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

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMVanGenuchten::WRMVanGenuchten(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};


/* ******************************************************************
* Relative permeability formula: input is capillary pressure pc.
* The original curve is regulized on interval (0, pc0) using the
* Hermite interpolant of order 3. Formulas (3.11)-(3.12).
****************************************************************** */
double WRMVanGenuchten::k_relative(double pc) {
  if (pc >= pc0_) {
    double se = pow(1.0 + pow(alpha_*pc, n_), -m_);
    if (function_ == FLOW_WRM_MUALEM) {
      return pow(se, l_) * pow(1.0 - pow(1.0 - pow(se, 1.0/m_), m_), 2.0);
    } else if (function_ == FLOW_WRM_BURDINE) {
      return se * se * (1.0 - pow(1.0 - pow(se, 1.0/m_), m_));
    } else if (function_ == FLOW_WRM_ONE) {
      return 1.;
    }
  } else if (pc <= 0.0) {
    return 1.0;
  } else {
    double pc_2 = pc * pc;
    double pc_3 = pc_2 * pc;
    return 1.0 + a_ * pc_2 + b_ * pc_3;
  }
  ASSERT(0);
  return 0.;
}


/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
double WRMVanGenuchten::d_k_relative(double pc) {
  if (pc >= pc0_) {
    double se = pow(1.0 + pow(alpha_*pc, n_), -m_);
    double dsdp = d_saturation(pc);

    double x = pow(se, 1.0 / m_);
    if (fabs(1.0 - x) < FLOW_WRM_TOLERANCE) return 0.0;

    double y = pow(1.0 - x, m_);
    double dkdse;
    if (function_ == FLOW_WRM_MUALEM)
      dkdse = (1.0 - y) * (l_ * (1.0 - y) + 2 * x * y / (1.0 - x)) * pow(se, l_ - 1.0);
    else if (function_ == FLOW_WRM_BURDINE)
      dkdse = (2 * (1.0 - y) + x / (1.0 - x)) * se;
    else if (function_ == FLOW_WRM_ONE)
      dkdse = 0.;

    double dk = dkdse * dsdp / (1 - sr_);
    ASSERT(std::abs(dk) < 1.e15);
    return dk;

  } else if (pc <= 0.0) {
    return 0.0;
  } else {
    return 2*a_*pc + 3*b_*pc*pc; 
  }
}


/* ******************************************************************
 * Saturation formula (3.5)-(3.6).
 ****************************************************************** */
double WRMVanGenuchten::saturation(double pc) {
  if (pc > 0.0) {
    return std::pow(1.0 + std::pow(alpha_*pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double WRMVanGenuchten::d_saturation(double pc) {
  if (pc > 0.0) {
    return -m_*n_ * std::pow(1.0 + std::pow(alpha_*pc, n_), -m_-1.0) * std::pow(alpha_*pc, n_-1) * alpha_ * (1.0 - sr_);
  } else {
    return 0.0;
  }
}

/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
double WRMVanGenuchten::capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return std::pow(se, -1.0/(m_*n_)) / alpha_;
  } else {
    return (std::pow(std::pow(se, -1.0/m_) - 1.0, 1/n_)) / alpha_;
  }
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double WRMVanGenuchten::d_capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return -1.0/(m_*n_*alpha_) * std::pow(se, -1.0/(m_*n_) - 1.)  / (1.0 - sr_);
  } else {
    return -1.0/(m_*n_*alpha_) * std::pow( std::pow(se, -1.0/m_) - 1.0, 1/n_-1.0)
        * std::pow(se, -1.0/m_ - 1.0) / (1.0 - sr_);
  }
}


void WRMVanGenuchten::InitializeFromPlist_() {
  std::string fname = plist_.get<std::string>("Krel function name", "Mualem");
  if (fname == std::string("Mualem")) {
    function_ = FLOW_WRM_MUALEM;
  } else if (fname == std::string("Burdine")) {
    function_ = FLOW_WRM_BURDINE;
  } else if (fname == std::string("one")) {
    function_ = FLOW_WRM_ONE;
  } else {
    ASSERT(0);
  }

  alpha_ = plist_.get<double>("van Genuchten alpha");
  sr_ = plist_.get<double>("residual saturation", 0.0);
  l_ = plist_.get<double>("Mualem exponent l", 0.5);
  pc0_ = plist_.get<double>("smoothing interval width", 0.0);

  if (plist_.isParameter("van Genuchten m")) {
    m_ = plist_.get<double>("van Genuchten m");
    if (function_ == FLOW_WRM_MUALEM || function_ == FLOW_WRM_ONE) {
      n_ = 1.0 / (1.0 - m_);
    } else {
      n_ = 2.0 / (1.0 - m_);
    }
  } else {
    n_ = plist_.get<double>("van Genuchten n");
    if (function_ == FLOW_WRM_MUALEM || function_ == FLOW_WRM_ONE) {
      m_ = 1.0 - 1.0/n_;
    } else {
      m_ = 1.0 - 2.0/n_;
    }
  }

  a_ = b_ = 0.;
  if (pc0_ > 0) {
    double k0 = k_relative(pc0_) - 1.0;
    double k0p = d_k_relative(pc0_);
    double pc0_2 = pc0_ * pc0_;
    double pc0_3 = pc0_2 * pc0_;

    a_ = (3 * k0 - k0p * pc0_) / pc0_2;
    b_ = (k0p * pc0_ - 2 * k0) / pc0_3;
  }

};

}  // namespace
}  // namespace
}  // namespace
