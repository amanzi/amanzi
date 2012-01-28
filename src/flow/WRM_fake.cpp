/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>

#include "WRM_fake.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRM_fake::WRM_fake(std::string region_)
{
  set_region(region_);
  alpha = 1.0;
  n = 2.0;
  m = 1.0;
}
 

/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRM_fake::k_relative(double pc)
{
  if (pc < 0.0) return 1.0 / (1.0 + pow(alpha*pc, n));
  else return 1.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov)                                     
****************************************************************** */
double WRM_fake::saturation(double pc)
{
  if (pc < 0.0) return pc;
  else return 1.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov).                                         
****************************************************************** */
double WRM_fake::d_saturation(double pc)
{
  if (pc < 0.0) return 1.0;
  else return 0.0;
}  


/* ******************************************************************
* erify (lipnikov@lanl.gov).                                       
****************************************************************** */
double WRM_fake::capillaryPressure(double sl)
{
  return 0.0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

