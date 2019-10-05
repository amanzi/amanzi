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
#include "errors.hh"
#include "WRM_BrooksCorey.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRM_BrooksCorey::WRM_BrooksCorey(const std::string region, double S_rw, double S_rn, double pd, double lambda)
{
  set_region(region);
  S_rw_ = S_rw;
  S_rn_ = S_rn;
  pd_ = pd;
  lambda_ = lambda;
}


/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
/*
double WRM_BrooksCorey::k_relative(double Sw, std::string phase_name)
{
  Errors::Message msg;
  double Swe = 0.0;
  Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  if (phase_name == "wetting") {
    return pow(Swe, 2.0);
  }
  else if (phase_name == "non wetting") {
    return pow(1.0 - Swe, 2.0);
  }
  else {
    msg << "Multiphase PK: phase_name \"" << phase_name.c_str() << "\" not recognized \n";
    Exceptions::amanzi_throw(msg);
  }
}
*/
double WRM_BrooksCorey::k_relative(double Sw, std::string phase_name)
{
  Errors::Message msg;
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  double Sne = 1.0 - Swe;
  if (phase_name == "wetting") {
    if (Swe < -1e-12) {
      return 0.0;
    } else if (Swe > 1.0) {
      return 1.0;
    } else {
      return pow(Swe,(2.0+3.0*lambda_)/lambda_);
    }
  }
  else if (phase_name == "non wetting") {
    if (Swe < -1e-12) {
      return 1.0;
    } else if (Swe > 1.0) {
      return 0.0;
    } else {
      return pow(Sne,2.0)*(1.0 - pow(Swe,(2+lambda_)/lambda_));
    }
  }
  else {
    msg << "Multiphase PK: phase_name \"" << phase_name.c_str() << "\" not recognized \n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Derivative of Relative permeability wrt wetting saturation formula.                                          
****************************************************************** */
/*
double WRM_BrooksCorey::dKdS(double Sw, std::string phase_name)
{
  Errors::Message msg;
  double Swe = 0.0;
  double factor = 1.0/(1.0 - S_rw_ - S_rn_);
  Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  if (phase_name == "wetting") {
    return 2.0 * Swe * factor;
  }
  else if (phase_name == "non wetting") {
    return - 2.0 * (1.0 - Swe) * factor;
  }
  else {
    msg << "Multiphase PK: phase_name \"" << phase_name.c_str() << "\" not recognized \n";
    Exceptions::amanzi_throw(msg);
  }
}
*/

double WRM_BrooksCorey::dKdS(double Sw, std::string phase_name)
{
  Errors::Message msg;
  double factor = 1.0/(1.0 - S_rw_ - S_rn_);
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  if (phase_name == "wetting") {
    return factor*(2.0 + 3.0*lambda_)/lambda_ * pow(Swe, (2.0 + 3.0*lambda_)/lambda_ - 1.0);
  }
  else if (phase_name == "non wetting") {
    return -(lambda_ + 2.0)/lambda_*factor*pow(1.0-Swe,2.0)*pow(Swe,(lambda_+2.0)/lambda_-1.0) - 
      2.0*factor*(1.0-Swe)*(1.0 - pow(Swe,(lambda_+2.0)/lambda_));
  }
  else {
    msg << "Multiphase PK: phase_name \"" << phase_name.c_str() << "\" not recognized \n";
    Exceptions::amanzi_throw(msg);
  }
}


/* *********************************************************************
* Capillary Pressure formula. Cut off capillary pressure for small Swe
***********************************************************************/
double WRM_BrooksCorey::capillaryPressure(double Sw)
{
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  if (Swe > 1e-12) {
    return pd_ * pow(Swe, -1.0/lambda_);
  } else {
    return pd_ * 20.0;
  }
}


/* ******************************************************************
* Derivative of capillary pressure formula.
****************************************************************** */
double WRM_BrooksCorey::dPc_dS(double Sw)
{
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  double factor = 1.0/(1.0 - S_rw_ - S_rn_);
  if (Swe > 1e-12) {
    return pd_ * factor * (-1.0/lambda_) * pow(Swe, -1.0/lambda_ - 1.0); 
  } else {
    return 1e5;
  }
}


/* ******************************************************************
* Return irreducible residual saturation of the phase.                                          
****************************************************************** */
double WRM_BrooksCorey::residualSaturation(std::string phase_name)
{ 
  if (phase_name == "wetting") {
    return S_rw_;
  } else if (phase_name == "non wetting") {
    return S_rn_;
  } else {
    Errors::Message msg;
    msg << "Unknown phase " << phase_name << ". Options are <wetting> and <non wetting>.\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi

