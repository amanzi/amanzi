/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>

#include "vanGenuchtenModel.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
vanGenuchtenModel::vanGenuchtenModel(
   std::string region_, double m_, double alpha_, double sr_, double atm_pressure_) :
   m(m_), alpha(alpha_), sr(sr_), atm_pressure(atm_pressure_)
{
  n = 1.0 / (1.0 - m);
  set_region(region_);
}
 

/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double vanGenuchtenModel::k_relative(double p)
{
  double pc = atm_pressure - p;  // capillary pressure
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
double vanGenuchtenModel::saturation(double p)
{
  double pc = atm_pressure - p;  // capillary pressure
  if (pc > 0.0) {
    return pow(1.0 + pow(alpha*pc, n), -m) * (1.0 - sr) + sr;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Derivative of the saturation formula w.r.t. capillary pressure.                                         
****************************************************************** */
double vanGenuchtenModel::d_saturation(double p)
{
  double pc = atm_pressure - p; // capillary pressure
  if (pc > 0.0) {
    return m*n * pow(1.0 + pow(alpha*pc, n), -m-1.0) * pow(alpha*pc, n-1) * alpha * (1.0 - sr);
  } else {
    return 0.0;
  }
}  


/* ******************************************************************
* Pressure as a function of saturation.                                       
****************************************************************** */
double vanGenuchtenModel::pressure(double sl)
{
  double se = (sl - sr) / (1.0 - sr);
  return atm_pressure - (pow(pow(se, -1.0/m) - 1.0, 1/n)) / alpha;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

