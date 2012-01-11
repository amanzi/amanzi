/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>

#include "WRM_analytic.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRM_analytic::WRM_analytic(std::string region_)
{
  set_region(region_);
  alpha = 1.0;
  atm_pressure = 0.0;
  n = 2.0;
  m = 1.0;
}
 

/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRM_analytic::k_relative(double p)
{
  if (p > 0.0) return 1.0 + pow(alpha*p, n);
  else return 1.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov)                                     
****************************************************************** */
double WRM_analytic::saturation(double p)
{
  if (p > 0.0) return 1.0 + pow(alpha*p, n);
  else return 1.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov).                                         
****************************************************************** */
double WRM_analytic::d_saturation(double p)
{
  if (p > 0.0) return 1.0;
  else return 0.0;
}  


/* ******************************************************************
* erify (lipnikov@lanl.gov).                                       
****************************************************************** */
double WRM_analytic::pressure(double sl)
{
  return 0.0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

