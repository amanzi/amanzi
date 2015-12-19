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
    return 1.0 / (1.0 + pow(alpha*pc, n));
  else
    return 1.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov)                                     
****************************************************************** */
double WRM_fake::saturation(double pc) const
{
  if (pc < 0.0)
    return 1.0; // 1.0 / (1.0 + pc * pc);
  else
    return 1.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov).                                         
****************************************************************** */
double WRM_fake::dSdPc(double pc) const
{
  if (pc < 0.0)
    return 0.0;  // pc / (1.0 + pc * pc);
  else
    return 0.0;
}


/* ******************************************************************
* Verify (lipnikov@lanl.gov).                                       
****************************************************************** */
double WRM_fake::capillaryPressure(double sl) const
{
  return 0.0;
}

}  // namespace Flow
}  // namespace Amanzi

