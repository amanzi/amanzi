/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <cmath>
#include <string>

#include "errors.hh"
#include "FlowDefs.hh"
#include "WRM_linear.hh"

namespace Amanzi {
namespace Flow {


/* ******************************************************************
* Constructors
****************************************************************** */
WRM_linear::WRM_linear(Teuchos::ParameterList& plist)
{
  pc0_ = plist.get<double>("pc0");
  sr_ = plist.get<double>("residual saturation liquid", FLOW_WRM_EXCEPTION);
}


WRM_linear::WRM_linear(double pc0) : pc0_(pc0) {}


/* ******************************************************************
* Relative permeability formula
****************************************************************** */
double
WRM_linear::k_relative(double pc) const
{
  return saturation(pc);
}


/* ******************************************************************
* Saturation formula
****************************************************************** */
double
WRM_linear::saturation(double pc) const
{
  if (pc < 0.0)
    return 1.0;
  else if (pc > pc0_)
    return sr_;
  else
    return sr_ + (1.0 - sr_) * (pc0_ - pc) / pc0_;
}


/* ******************************************************************
* Pressure as a function of saturation.
****************************************************************** */
double
WRM_linear::capillaryPressure(double s) const
{
  return pc0_ * (1.0 - s) / (1.0 - sr_);
}


/* ******************************************************************
* Derivatives
****************************************************************** */
double
WRM_linear::dSdPc(double pc) const
{
  if (pc < 0.0)
    return 0.0;
  else if (pc > pc0_)
    return 0.0;
  return -(1.0 - sr_) / pc0_;
}


double
WRM_linear::dKdPc(double pc) const
{
  return dSdPc(pc);
}

} // namespace Flow
} // namespace Amanzi
