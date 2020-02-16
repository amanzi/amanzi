/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#include <cmath>
#include <string>
#include <iostream>

// Amanzi
#include "errors.hh"

// Amanzi::Multiphase
#include "MultiphaseDefs.hh"
#include "WRMmp_vanGenuchten.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRMmp_vanGenuchten::WRMmp_vanGenuchten(Teuchos::ParameterList& plist)
{
  double srl = plist.get<double>("residual saturation liquid", MULTIPHASE_WRM_EXCEPTION);
  double srg = plist.get<double>("residual saturation gas", MULTIPHASE_WRM_EXCEPTION);
  double n = plist.get<double>("van Genuchten n", MULTIPHASE_WRM_EXCEPTION);
  double Pr = plist.get<double>("van Genuchten entry pressure", MULTIPHASE_WRM_EXCEPTION);
  double reg_kr = plist.get<double>("regularization interval kr", MULTIPHASE_WRM_REGULARIZATION_INTERVAL);
  double reg_pc = plist.get<double>("regularization interval pc", MULTIPHASE_WRM_REGULARIZATION_INTERVAL);

  Init_(srl, srg, n, Pr, reg_kr, reg_pc);
}

void WRMmp_vanGenuchten::Init_(double srl, double srg, double n, double Pr,
                               double reg_kr, double reg_pc)
{
  srl_ = srl;
  srg_ = srg;
  n_ = n;
  m_ = 1.0 - 1.0 / n_;
  Pr_ = Pr;
  reg_kl_ = reg_kr;
  reg_pc_ = reg_pc;

  if (reg_kl_ > 0.0) {
    double se0(1.0 - reg_kl_), se1(1.0);
    spline_kl_.Setup(se0, k_relative_liquid_(se0), dKdS_liquid_(se0),
                     se1, 1.0, 0.0);

    if (spline_kl_.grad().NormInf() > MULTIPHASE_WRM_REGULARIZATION_MAX_GRADIENT) {
      Errors::Message msg;
      msg << "Gradient of regularization spline exceeds threshold " << MULTIPHASE_WRM_REGULARIZATION_MAX_GRADIENT;
      Exceptions::amanzi_throw(msg);
    }
  }

  if (reg_pc_ > 0.0) {
    double se0(reg_pc_), se1(1.0 - reg_pc_);
    spline_pc_.Setup(se0, capillaryPressure_(se0), dPc_dS_(se0),
                     se1, capillaryPressure_(se1), dPc_dS_(se1));
  }
}


/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRMmp_vanGenuchten::k_relative(double sl, const std::string& phase)
{
  double sle = (sl - srl_) / (1.0 - srl_ - srg_);
  if (phase == "liquid") {
    if (sle <= 0.0) {
      return 0.0;
    } else if (sle >= 1.0) {
      return 1.0;
    } else if (sle > 1.0 - reg_kl_) {
      return spline_kl_.Value(sle);
    } else {
      return k_relative_liquid_(sle);
    }
  }
  else if (phase == "gas") {
    if (sle <= 0.0) {
      return 1.0;
    } else if (sle >= 1.0) {
      return 0.0;
    } else {
      return k_relative_gas_(sle);
    }
  }
}


/* ******************************************************************
* Relative permeability wrt to effective saturation.
****************************************************************** */
double WRMmp_vanGenuchten::k_relative_liquid_(double sle) {
  return pow(sle, 0.5) * pow(1.0 - pow(1.0 - pow(sle, 1.0 / m_), m_), 2.0);
}

double WRMmp_vanGenuchten::k_relative_gas_(double sle) {
  return pow(1.0 - sle, 0.5) * pow(1.0 - pow(sle, 1.0 / m_), 2.0 * m_);
}


/* ******************************************************************
* Derivative of relative permeability wrt liquid saturation
****************************************************************** */
double WRMmp_vanGenuchten::dKdS(double sl, const std::string& phase)
{
  double sle = (sl - srl_) / (1.0 - srl_ - srg_);
  if (phase == "liquid") {
    if (sle <= 0.0) {
      return 0.0;
    } else if (sle >= 1.0) {
      return 0.0;
    } else if (sle > 1.0 - reg_kl_) {
      return spline_kl_.GradientValue(sle);
    } else {
      return dKdS_liquid_(sle);
    }
  }
  else if (phase == "gas") {
    if (sle <= 0.0) {
      return 0.0;
    } else if (sle >= 1.0) {
      return 0.0;
    } else {
      return dKdS_gas_(sle);
    }
  }
}


/* ******************************************************************
* Derivative of relative permeability wrt effective saturation
****************************************************************** */
double WRMmp_vanGenuchten::dKdS_liquid_(double sle) {
  double factor = 1.0 / (1.0 - srl_ - srg_);
  double a1 = std::pow(sle, 1.0 / m_);
  double a2 = 1.0 - a1;
  double a3 = std::pow(a2, m_ - 1.0);
  double a4 = 1.0 - a3 * a2;
  double a5 = a4 * (a4 / 2 + 2 * a3 * a1);

  return factor * a5 / std::pow(sle, 0.5);
}

double WRMmp_vanGenuchten::dKdS_gas_(double sle) {
  double factor = 1.0 / (1.0 - srl_ - srg_);
  double a1 = std::pow(sle, 1.0 / m_);
  double a2 = 1.0 - a1;
  double a3 = std::pow(a2, 2 * m_ - 1.0);
  double a4 = std::pow(a2, 2 * m_);
  double b1 = 2 * std::pow(1.0 - sle, 0.5);

  return -factor * (a4 / b1 + a3 * b1 * std::pow(sle, 1.0/m_ - 1.0));
}


/* ******************************************************************
* Capillary Pressure formula. 
****************************************************************** */
double WRMmp_vanGenuchten::capillaryPressure(double sl)
{
  double sle = (sl - srl_) / (1.0 - srl_ - srg_);
  if (sle <= reg_pc_) {
    return spline_pc_.Value(sle);
  } else if (sle >= 1.0 - reg_pc_) {
    return spline_pc_.Value(sle);
  } else {
    return capillaryPressure_(sle);
  }
}


/* ******************************************************************
* Derivative of capillary pressure formula. 
****************************************************************** */
double WRMmp_vanGenuchten::dPc_dS(double sl)
{
  double factor = 1.0 / (1.0 - srl_ - srg_);
  double sle = (sl - srl_) / (1.0 - srl_ - srg_);
  if (sle <= reg_pc_) {
    return factor * spline_pc_.GradientValue(sle);
  } else if (sle >= 1.0 - reg_pc_) {
    return factor * spline_pc_.GradientValue(sle);
  } else {
    return dPc_dS_(sle);
  }
}


/* ******************************************************************
* Return irreducible residual saturation of the phase.                                          
****************************************************************** */
double WRMmp_vanGenuchten::capillaryPressure_(double sle) {
  return Pr_ * pow(pow(sle, -1.0 / m_) - 1.0, 1.0 / n_);
}


double WRMmp_vanGenuchten::dPc_dS_(double sle) {
  return -Pr_ / (m_ * n_) * pow(pow(sle, -1.0/m_) - 1.0, 1.0/n_ - 1.0) * pow(sle, -1.0/m_ - 1.0) / (1.0 - srl_ - srg_);
}

}  // namespace Multiphase
}  // namespace Amanzi

