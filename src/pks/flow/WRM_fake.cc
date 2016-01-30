/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

#include "WRM_fake.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Setup fundamental parameters for this model.                                            
****************************************************************** */
WRM_fake::WRM_fake(Teuchos::ParameterList& plist)
{
  set_region("dummy");
  alpha = 1.0;
  n = 2.0;
  m = 1.0;
}


/* ******************************************************************
* Relative permeability formula.                                          
****************************************************************** */
double WRM_fake::k_relative(double pc) const
{
  if (pc < 0.0)
    return 1.0 / (1.0 + pc * pc);
  else
    return 1.0;
}


/* ******************************************************************
* The analytic solution was designed for stationary PDE.
****************************************************************** */
double WRM_fake::saturation(double pc) const {
  return -pc;
}


/* ******************************************************************
* Derivative of rel_perm  w.r.t. capillary pressure.                                     
****************************************************************** */
double WRM_fake::dKdPc(double pc) const
{
  if (pc < 0.0) {
    double tmp = 1.0 + pc * pc;
    return -2 * pc / (tmp * tmp);
  } else {
    return 0.0;
  }
}


/* ******************************************************************
* The analytic solution was designed for stationary PDE.
****************************************************************** */
double WRM_fake::dSdPc(double pc) const {
  return -1.0; 
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov).                                       
****************************************************************** */
double WRM_fake::capillaryPressure(double sl) const {
  return -sl;
}

}  // namespace Flow
}  // namespace Amanzi

