/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

#include "WRM_vanGenuchten.hpp"

namespace Amanzi {
namespace AmanziFlow {

const int FLOW_WRM_MUALEM = 1;
const int FLOW_WRM_BURDINE = 2;

/* ******************************************************************
* Setup fundamental parameters for this model.
* Default value of the regularization interval is pc0 = 0.                                           
****************************************************************** */
WRM_vanGenuchten::WRM_vanGenuchten(
    std::string region, double m, double l, double alpha, 
    double sr, std::string krel_function, double pc0)
    : m_(m), l_(l), alpha_(alpha), sr_(sr), pc0_(pc0)
{
  set_region(region);

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
double WRM_vanGenuchten::k_relative(double pc)
{
  if (pc >= pc0_) {
    double se = pow(1.0 + pow(alpha_*pc, n_), -m_);
    if (function_ == FLOW_WRM_MUALEM) {
      return pow(se, l_) * pow(1.0 - pow(1.0 - pow(se, 1.0/m_), m_), 2.0);
    } else {
      return se * se * (1.0 - pow(1.0 - pow(se, 1.0/m_), m_));     
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
double WRM_vanGenuchten::saturation(double pc)
{
  if (pc > 0.0) {
    return pow(1.0 + pow(alpha_*pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Derivative of the saturation formula w.r.t. capillary pressure.
* Warning: remember that dSdP = -dSdPc.                                        
****************************************************************** */
double WRM_vanGenuchten::dSdPc(double pc)
{
  if (pc > 0.0) {
    double alpha_pc = alpha_*pc;
    double x = pow(alpha_pc, n_-1.0);
    double y = x * alpha_pc;
    return pow(1.0 + y, -m_-1.0) * x * factor_dSdPc_;
  } else {
    return 0.0;
  }
}


/* ******************************************************************
* Pressure as a function of saturation.                                       
****************************************************************** */
double WRM_vanGenuchten::capillaryPressure(double s)
{
  double se = (s - sr_) / (1.0 - sr_);
  return (pow(pow(se, -1.0/m_) - 1.0, 1/n_)) / alpha_;
}


/* ******************************************************************
* Derivative of the original relative permeability w.r.t. capillary pressure.                                     
****************************************************************** */
double WRM_vanGenuchten::dKdPc(double pc)
{
  if (pc > 0.0) {
    double se = pow(1.0 + pow(alpha_*pc, n_), -m_);
    double dsdp = dSdPc(pc);

    double x = pow(se, 1.0 / m_);
    double y = pow(1.0 - x, m_);
    double dkdse;
    if (function_ == FLOW_WRM_MUALEM)
      dkdse = (1.0 - y) * (l_ * (1.0 - y) + 2 * x * y / (1.0 - x)) * pow(se, l_ - 1.0);
    else
      dkdse = (2 * (1.0 - y) + x / (1.0 - x)) * se; 

    return dkdse * dsdp / (1 - sr_);

  } else {
    return 0.0;
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

