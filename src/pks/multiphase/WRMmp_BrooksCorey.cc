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
#include "MultiphaseDefs.hh"
#include "WRMmp_BrooksCorey.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRMmp_BrooksCorey::WRMmp_BrooksCorey(Teuchos::ParameterList& plist)
{
  double S_rw = plist.get<double>("residual saturation liquid", MULTIPHASE_WRM_EXCEPTION);
  double S_rn = plist.get<double>("residual saturation gas", MULTIPHASE_WRM_EXCEPTION);
  double pd = plist.get<double>("entry pressure", MULTIPHASE_WRM_EXCEPTION);
  double lambda = plist.get<double>("brooks corey lambda", MULTIPHASE_WRM_EXCEPTION);

  Init_(S_rw, S_rn, pd, lambda);
}

void WRMmp_BrooksCorey::Init_(double S_rw, double S_rn, double pd, double lambda)
{
  S_rw_ = S_rw;
  S_rn_ = S_rn;
  pd_ = pd;
  lambda_ = lambda;
}


double WRMmp_BrooksCorey::k_relative(double Sw, int phase)
{
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  double Sne = 1.0 - Swe;
  if (phase == MULTIPHASE_PHASE_LIQUID) {
    if (Swe < -1e-12) {
      return 0.0;
    } else if (Swe > 1.0) {
      return 1.0;
    } else {
      return pow(Swe,(2.0+3.0*lambda_)/lambda_);
    }
  }
  else if (phase == MULTIPHASE_PHASE_GAS) {
    if (Swe < -1e-12) {
      return 1.0;
    } else if (Swe > 1.0) {
      return 0.0;
    } else {
      return pow(Sne,2.0)*(1.0 - pow(Swe,(2+lambda_)/lambda_));
    }
  }
}


/* *********************************************************************
* TBW
***********************************************************************/
double WRMmp_BrooksCorey::dKdS(double Sw, int phase)
{
  double factor = 1.0/(1.0 - S_rw_ - S_rn_);
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  if (phase == MULTIPHASE_PHASE_LIQUID) {
    return factor*(2.0 + 3.0*lambda_)/lambda_ * pow(Swe, (2.0 + 3.0*lambda_)/lambda_ - 1.0);
  }
  else if (phase == MULTIPHASE_PHASE_GAS) {
    return -(lambda_ + 2.0)/lambda_*factor*pow(1.0-Swe,2.0)*pow(Swe,(lambda_+2.0)/lambda_-1.0) - 
      2.0*factor*(1.0-Swe)*(1.0 - pow(Swe,(lambda_+2.0)/lambda_));
  }
}


/* *********************************************************************
* Capillary Pressure formula. Cut off capillary pressure for small Swe
***********************************************************************/
double WRMmp_BrooksCorey::capillaryPressure(double Sw)
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
double WRMmp_BrooksCorey::dPc_dS(double Sw)
{
  double Swe = (Sw - S_rw_)/(1.0 - S_rw_ - S_rn_);
  double factor = 1.0/(1.0 - S_rw_ - S_rn_);
  if (Swe > 1e-12) {
    return pd_ * factor * (-1.0/lambda_) * pow(Swe, -1.0/lambda_ - 1.0); 
  } else {
    return 1e5;
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

