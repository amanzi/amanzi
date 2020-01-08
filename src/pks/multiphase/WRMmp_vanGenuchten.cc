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

#include "errors.hh"
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

  Init_(srl, srg, n, Pr);
}

void WRMmp_vanGenuchten::Init_(double srl, double srg, double n, double Pr)
{
  srl_ = srl;
  srg_ = srg;
  n_ = n;
  m_ = 1.0 - 1.0/n_;
  Pr_ = Pr;
  eps_ = 1.0e-3;
}


/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRMmp_vanGenuchten::k_relative(double Sw, std::string phase_name)
{
  Errors::Message msg;
  double Swe = (Sw - srl_) / (1.0 - srl_ - srg_);
  if (phase_name == "liquid") {
    if (Swe < 1.0e-09) {
      return 0.0;
    } else if(Swe - 1.0 > -1.0e-09) {
      return 1.0;
    } else {
      double tmp = pow(Swe, 0.5) * pow(1.0 - pow(1.0 - pow(Swe, 1.0/m_), m_), 2.0);
      return tmp; 
    }
  }
  else if (phase_name == "gas") {
    if (Swe < 1.0e-09) {
      return 1.0;
    } else if(Swe - 1.0 > -1.0e-09) {
      return 0.0;
    } else {
      return pow(1.0 - Swe, 0.5) * pow(1.0 - pow(Swe, 1.0/m_), 2.0 * m_);
    }
  }
  else {
    msg << "Multiphase PK: phase_name \"" << phase_name.c_str() << "\" not recognized \n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Derivative of relative permeability wrt liquid saturation
****************************************************************** */
double WRMmp_vanGenuchten::dKdS(double Sw, std::string phase_name)
{
  Errors::Message msg;
  double Swe = 0.0;
  double factor = 1.0/(1.0 - srl_ - srg_);
  Swe = (Sw - srl_)/(1.0 - srl_ - srg_);
  if (phase_name == "liquid") {
    if (Swe < 1.0e-09) {
      return 0.0;
    } else if(Swe - 1.0 > -1.0e-09) {
      return 0.5 * factor;
    } else {
      return factor*0.5*pow(Swe,-0.5)*pow(1.0 - pow(1.0 - pow(Swe,1.0/m_),m_),2.0) + 
        2.0*(1.0 - pow(1.0 - pow(Swe, 1.0/m_),m_))*pow(1.0 - pow(Swe, 1.0/m_), m_ - 1.0)*pow(Swe, 1.0/m_ - 0.5)*factor;
    }
  }
  else if (phase_name == "gas") {
    if (Swe < 1.0e-09) {
      return -0.5 * factor;
    } else if(Swe - 1.0 > -1.0e-09) {
      return 0.0;
    } else {
      return -factor*0.5*pow(1.0 - Swe,-0.5)*pow(1.0 - pow(Swe, 1.0/m_), 2.0*m_) - 
      factor*pow(1.0 - Swe, 0.5)*2.0*pow(1.0 - pow(Swe, 1.0/m_), 2.0*m_ - 1.0)*pow(Swe, 1.0/m_ - 1.0);
    }
  }
  else {
    msg << "Multiphase PK: phase_name \"" << phase_name.c_str() << "\" not recognized \n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Capillary Pressure formula. 
****************************************************************** */
double WRMmp_vanGenuchten::capillaryPressure(double Sw)
{
  double Sn = 1.0 - Sw;
  if (Sn - srg_ < 1e-15) {
    return mod_VGM(srg_) + deriv_mod_VGM(srg_) * (Sn - srg_);
  } else if (Sn + srl_ - 1.0 > -1e-15) {
    return mod_VGM(1.0 - srl_) + deriv_mod_VGM(1.0 - srl_) * (Sn - 1.0 + srl_);
  } else {
    return mod_VGM(Sn);
  }
}


/* ******************************************************************
* Derivative of capillary pressure formula. 
****************************************************************** */
double WRMmp_vanGenuchten::dPc_dS(double Sw)
{
  double Sn = 1.0 - Sw;
  if (Sn - srg_ < 1e-15) {
    return -deriv_mod_VGM(srg_);
  } else if (Sn + srl_ - 1.0 > -1e-15) {
    return -deriv_mod_VGM(1.0 - srl_);
  } else {
    return -deriv_mod_VGM(Sn); // wrt Sw
  }
}



/* ******************************************************************
* Return irreducible residual saturation of the phase.                                          
****************************************************************** */
double WRMmp_vanGenuchten::VGM(double Sn) {
  double Sw = 1.0 - Sn;
  double Swe = (Sw - srl_)/(1.0 - srl_ - srg_);
  return Pr_ * pow(pow(Swe, -1.0/m_) - 1.0, 1.0/n_);
}

double WRMmp_vanGenuchten::mod_VGM(double Sn) {
  double s_mod = srg_ + (1.0 - eps_) * (Sn - srg_) + eps_/2.0 * (1.0 - srl_ - srg_);
  return VGM(s_mod) - VGM(srg_ + eps_/2 * (1.0 - srl_ - srg_));
}

double WRMmp_vanGenuchten::deriv_VGM(double Sn) { // wrt Sn
  double Sw = 1.0 - Sn;
  double Swe = (Sw - srl_)/(1.0 - srl_ - srg_);
  return Pr_ / (m_ * n_) * pow(pow(Swe, -1.0/m_) - 1.0, 1.0/n_ - 1.0) * pow(Swe, -1.0/m_ - 1.0) / (1.0 - srl_ - srg_);
}

double WRMmp_vanGenuchten::deriv_mod_VGM(double Sn) { // wrt Sn
  double s_mod = srg_ + (1.0 - eps_) * (Sn - srg_) + eps_/2.0 * (1.0 - srl_ - srg_);
  return (1.0 - eps_) * deriv_VGM(s_mod);
}

}  // namespace Multiphase
}  // namespace Amanzi

