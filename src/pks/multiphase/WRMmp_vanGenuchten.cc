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
#include "VectorObjects.hh"

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
  double reg = plist.get<double>("regularization interval", MULTIPHASE_WRM_REGULARIZATION_INTERVAL);

  Init_(srl, srg, n, Pr, reg);
}

void WRMmp_vanGenuchten::Init_(double srl, double srg, double n, double Pr, double reg)
{
  srl_ = srl;
  srg_ = srg;
  n_ = n;
  m_ = 1.0 - 1.0 / n_;
  Pr_ = Pr;
  eps_ = 1.0e-3;

  reg_kl_ = reg;
  if (reg_kl_ > 0.0) {
    double se0(1.0 - reg_kl_), se1(1.0);
    spline_kl_.Setup(se0, k_relative_liquid_(se0), dKdS_liquid_(se0),
                     se1, 1.0, 0.0);
    grad_spline_kl_ = (WhetStone::Gradient(spline_kl_.poly()))[0];

    if (grad_spline_kl_.NormInf() > MULTIPHASE_WRM_REGULARIZATION_MAX_GRADIENT) {
      Errors::Message msg;
      msg << "Gradient of regularization spline exceeds threshold " << MULTIPHASE_WRM_REGULARIZATION_MAX_GRADIENT;
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRMmp_vanGenuchten::k_relative(double sl, const std::string& phase)
{
  if (phase == "liquid") {
    double sle = (sl - srl_) / (1.0 - srl_);
    if (sle < 0.0) {
      return 0.0;
    } else if (sle > 1.0 - reg_kl_) {
      return spline_kl_.Value(sle);
    } else {
      return k_relative_liquid_(sle);
    }
  }
  else if (phase == "gas") {
    double sle = (sl - srl_) / (1.0 - srl_ - srg_);
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
  if (phase == "liquid") {
    double sle = (sl - srl_) / (1.0 - srl_);
    if (sle < 0.0) {
      return 0.0;
    } else if (sle > 1.0 - reg_kl_) {
      AmanziGeometry::Point x(1);
      x[0] = sle;
      return grad_spline_kl_.Value(x);
    } else {
      return dKdS_liquid_(sle);
    }
  }
  else if (phase == "gas") {
    double sle = (sl - srl_) / (1.0 - srl_ - srg_);
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
  double sg = 1.0 - sl;
  if (sg - srg_ < 1e-15) {
    return mod_VGM(srg_) + deriv_mod_VGM(srg_) * (sg - srg_);
  } else if (sg + srl_ - 1.0 > -1e-15) {
    return mod_VGM(1.0 - srl_) + deriv_mod_VGM(1.0 - srl_) * (sg - 1.0 + srl_);
  } else {
    return mod_VGM(sg);
  }
}


/* ******************************************************************
* Derivative of capillary pressure formula. 
****************************************************************** */
double WRMmp_vanGenuchten::dPc_dS(double sl)
{
  double sg = 1.0 - sl;
  if (sg - srg_ < 1e-15) {
    return -deriv_mod_VGM(srg_);
  } else if (sg + srl_ - 1.0 > -1e-15) {
    return -deriv_mod_VGM(1.0 - srl_);
  } else {
    return -deriv_mod_VGM(sg);  // negative wrt sl
  }
}


/* ******************************************************************
* Return irreducible residual saturation of the phase.                                          
****************************************************************** */
double WRMmp_vanGenuchten::VGM(double sg) {
  double sl = 1.0 - sg;
  double sle = (sl - srl_) / (1.0 - srl_ - srg_);
  return Pr_ * pow(pow(sle, -1.0 / m_) - 1.0, 1.0 / n_);
}


double WRMmp_vanGenuchten::mod_VGM(double sg) {
  double s_mod = srg_ + (1.0 - eps_) * (sg - srg_) + eps_/2.0 * (1.0 - srl_ - srg_);
  return VGM(s_mod) - VGM(srg_ + eps_/2 * (1.0 - srl_ - srg_));
}


double WRMmp_vanGenuchten::deriv_VGM(double sg) { // wrt sg
  double sl = 1.0 - sg;
  double sle = (sl - srl_)/(1.0 - srl_ - srg_);
  return Pr_ / (m_ * n_) * pow(pow(sle, -1.0/m_) - 1.0, 1.0/n_ - 1.0) * pow(sle, -1.0/m_ - 1.0) / (1.0 - srl_ - srg_);
}


double WRMmp_vanGenuchten::deriv_mod_VGM(double sg) { // wrt sg
  double s_mod = srg_ + (1.0 - eps_) * (sg - srg_) + eps_/2.0 * (1.0 - srl_ - srg_);
  return (1.0 - eps_) * deriv_VGM(s_mod);
}

}  // namespace Multiphase
}  // namespace Amanzi

