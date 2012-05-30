/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

#include "Flow_PK.hpp"
#include "WRM_BrooksCorey.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Setup fundamental parameters for this model.
* Default value of the regularization interval is pc0 = 0.                                           
****************************************************************** */
WRM_BrooksCorey::WRM_BrooksCorey(
    std::string region, double lambda, double alpha, double sr, std::string krel_function, double pc0)
    : lambda_(lambda), alpha_(alpha), sr_(sr), pc0_(pc0)
{
  set_region(region);
  if (krel_function == "Mualem") {
    factor_ = -2.0 - 2.5 * lambda_;
  } else if (krel_function == "Burdine") {
    factor_ = -2.0 - 3.0 * lambda_;
  }
}


/* ******************************************************************
* Relative permeability formula: input is capillary pressure pc.
* The original curve is regulized on interval (0, pc0) using the 
* Hermite interpolant of order 3.       
****************************************************************** */
double WRM_BrooksCorey::k_relative(double pc)
{
  if (pc > 0.0) {
    return pow(alpha_ * pc, -factor_);
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Saturation formula (3.5)-(3.8).                                         
****************************************************************** */
double WRM_BrooksCorey::saturation(double pc)
{
  if (pc > 0.0) {
    return pow(alpha_ * pc, -lambda_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Derivative of the saturation formula w.r.t. capillary pressure.
* Warning: remember that dSdP = -dSdPc.                                        
****************************************************************** */
double WRM_BrooksCorey::dSdPc(double pc)
{
  if (pc > 0.0) {
    return pow(alpha_ * pc, -lambda_ - 1.0) * (1.0 - sr_) * alpha_;
  } else {
    return 0.0;
  }
}


/* ******************************************************************
* Pressure as a function of saturation, formula (3.9).                                       
****************************************************************** */
double WRM_BrooksCorey::capillaryPressure(double s)
{
  double se = (s - sr_) / (1.0 - sr_);
  return pow(se, -1.0/alpha_) / alpha_;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

