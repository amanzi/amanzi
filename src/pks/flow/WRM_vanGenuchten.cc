/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
      Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <cmath>
#include <string>

#include "errors.hh"
#include "FlowDefs.hh"
#include "WRM_vanGenuchten.hh"

namespace Amanzi {
namespace Flow {

const int FLOW_WRM_MUALEM = 1;
const int FLOW_WRM_BURDINE = 2;
const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
* Setup fundamental parameters for this model.
* Default value of the regularization interval is pc0 = 0.
****************************************************************** */
WRM_vanGenuchten::WRM_vanGenuchten(Teuchos::ParameterList& plist)
{
  double m = plist.get<double>("van Genuchten m", FLOW_WRM_EXCEPTION);
  double alpha = plist.get<double>("van Genuchten alpha", FLOW_WRM_EXCEPTION);
  double l = plist.get<double>("van Genuchten l", FLOW_WRM_VANGENUCHTEN_L);
  double sr = plist.get<double>("residual saturation liquid", FLOW_WRM_EXCEPTION);
  double pc0 = plist.get<double>("regularization interval", FLOW_WRM_REGULARIZATION_INTERVAL);
  std::string krel_function = plist.get<std::string>("relative permeability model", "Mualem");

  Init_(m, l, alpha, sr, krel_function, pc0);
}


WRM_vanGenuchten::WRM_vanGenuchten(double m,
                                   double l,
                                   double alpha,
                                   double sr,
                                   std::string& krel_function,
                                   double pc0)
{
  Init_(m, l, alpha, sr, krel_function, pc0);
}


/* ******************************************************************
* Setup fundamental parameters for this model.
* Default value of the regularization interval is pc0 = 0.
****************************************************************** */
void
WRM_vanGenuchten::Init_(double m,
                        double l,
                        double alpha,
                        double sr,
                        std::string& krel_function,
                        double pc0)
{
  m_ = m;
  l_ = l;
  alpha_ = alpha;
  sr_ = sr;
  pc0_ = pc0;
  tol_ = FLOW_WRM_TOLERANCE;

  Errors::Message msg;
  if (m_ < 0.0 || alpha_ < 0.0 || sr_ < 0.0 || pc0_ < 0.0) {
    msg << "vanGenuchten: negative parameter in a water retention model.";
    Exceptions::amanzi_throw(msg);
  }
  if (sr_ > 1.0) {
    msg << "vanGenuchten: residual saturation is greater than 1.";
    Exceptions::amanzi_throw(msg);
  }

  if (krel_function == "Mualem") {
    n_ = 1.0 / (1.0 - m_);
    function_ = FLOW_WRM_MUALEM;
  } else {
    n_ = 2.0 / (1.0 - m_);
    function_ = FLOW_WRM_BURDINE;
  }
  factor_dSdPc_ = -m_ * n_ * alpha_ * (1.0 - sr_);
  a_ = b_ = 0;

  if (pc0 > 0) {
    double k0 = k_relative(pc0) - 1.0;
    double k0p = dKdPc(pc0);
    double pc0_2 = pc0 * pc0;
    double pc0_3 = pc0_2 * pc0;

    a_ = (3 * k0 - k0p * pc0) / pc0_2;
    b_ = (k0p * pc0 - 2 * k0) / pc0_3;
  }
}


/* ******************************************************************
* Relative permeability formula: input is capillary pressure pc.
* The original curve is regulized on interval (0, pc0) using the
* Hermite interpolant of order 3. Formulas (3.11)-(3.12).
****************************************************************** */
double
WRM_vanGenuchten::k_relative(double pc) const
{
  if (pc >= pc0_) {
    double se = pow(1.0 + pow(alpha_ * pc, n_), -m_);
    if (function_ == FLOW_WRM_MUALEM) {
      double tmp = 1.0 - pow(1.0 - pow(se, 1.0 / m_), m_);
      return pow(se, l_) * tmp * tmp;
    } else {
      return se * se * (1.0 - pow(1.0 - pow(se, 1.0 / m_), m_));
    }
  } else if (pc <= 0.0) {
    return 1.0;
  } else {
    double pc_2 = pc * pc;
    double pc_3 = pc_2 * pc;
    return 1.0 + a_ * pc_2 + b_ * pc_3;
  }
}


/* ******************************************************************
* Saturation formula (3.5)-(3.6).
****************************************************************** */
double
WRM_vanGenuchten::saturation(double pc) const
{
  if (pc > 0.0) {
    return pow(1.0 + pow(alpha_ * pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Derivative of the saturation formula w.r.t. capillary pressure.
* Warning: remember that dSdP = -dSdPc.
****************************************************************** */
double
WRM_vanGenuchten::dSdPc(double pc) const
{
  if (pc > 0.0) {
    double alpha_pc = alpha_ * pc;
    double x = pow(alpha_pc, n_ - 1.0);
    double y = x * alpha_pc;
    return pow(1.0 + y, -m_ - 1.0) * x * factor_dSdPc_;
  } else {
    return 0.0;
  }
}


/* ******************************************************************
* Pressure as a function of saturation.
****************************************************************** */
double
WRM_vanGenuchten::capillaryPressure(double s) const
{
  double se = (s - sr_) / (1.0 - sr_);
  return (pow(pow(se, -1.0 / m_) - 1.0, 1.0 / n_)) / alpha_;
}


/* ******************************************************************
* Derivative of the original relative permeability w.r.t. capillary pressure.
****************************************************************** */
double
WRM_vanGenuchten::dKdPc(double pc) const
{
  if (pc >= pc0_) {
    double se = pow(1.0 + pow(alpha_ * pc, n_), -m_);
    double dsdp = dSdPc(pc);

    double x = pow(se, 1.0 / m_);
    if (fabs(1.0 - x) < tol_) return 0.0;

    double y = pow(1.0 - x, m_);
    double dkdse;
    if (function_ == FLOW_WRM_MUALEM)
      dkdse = (1.0 - y) * (l_ * (1.0 - y) + 2 * x * y / (1.0 - x)) * pow(se, l_ - 1.0);
    else
      dkdse = (2 * (1.0 - y) + x / (1.0 - x)) * se;

    return dkdse * dsdp / (1.0 - sr_);

  } else if (pc <= 0.0) {
    return 0.0;

  } else {
    return 2 * a_ * pc + 3 * b_ * pc * pc;
  }
}

} // namespace Flow
} // namespace Amanzi
