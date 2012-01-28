/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>

#include "WRM_vanGenuchten.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRM_vanGenuchten::WRM_vanGenuchten(
   std::string region_, double m_, double alpha_, double sr_) :
   m(m_), alpha(alpha_), sr(sr_)
{
  n = 1.0 / (1.0 - m);
  set_region(region_);
}
 

/* ******************************************************************
* Relative permeability formula: input is capillary pressure.                                        
****************************************************************** */
double WRM_vanGenuchten::k_relative(double pc)
{
  if (pc > 0.0) {
    double se = pow(1.0 + pow(alpha*pc, n), -m);
    return sqrt(se) * pow(1.0 - pow(1.0 - pow(se, 1.0/m), m), 2.0);
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Saturation formula (3.5)-(3.6).                                         
****************************************************************** */
double WRM_vanGenuchten::saturation(double pc)
{
  if (pc > 0.0) {
    return pow(1.0 + pow(alpha*pc, n), -m) * (1.0 - sr) + sr;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Derivative of the saturation formula w.r.t. capillary pressure.                                         
****************************************************************** */
double WRM_vanGenuchten::d_saturation(double pc)
{
  if (pc > 0.0) {
    return m*n * pow(1.0 + pow(alpha*pc, n), -m-1.0) * pow(alpha*pc, n-1) * alpha * (1.0 - sr);
  } else {
    return 0.0;
  }
}  


/* ******************************************************************
* Pressure as a function of saturation.                                       
****************************************************************** */
double WRM_vanGenuchten::capillaryPressure(double sl)
{
  double se = (sl - sr) / (1.0 - sr);
  return (pow(pow(se, -1.0/m) - 1.0, 1/n)) / alpha;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

