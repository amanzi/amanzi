/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Milliton-Quirk constitutive model for tortuosity.
*/

#include <cmath>

#include "COM_MillingtonQuirk.hh"

namespace Amanzi {
namespace AmanziEOS {

double
COM_MillingtonQuirk::Tortuosity(double phi, double s)
{
  return pow(phi, a_) * pow(1.0 - s, b_);
}


double
COM_MillingtonQuirk::DTortuosityDphi(double phi, double s)
{
  return a_ * pow(phi, a_ - 1.0) * pow(1.0 - s, b_);
}


double
COM_MillingtonQuirk::DTortuosityDs(double phi, double s)
{
  return -b_ * pow(phi, a_) * pow(1.0 - s, b_ - 1);
}

} // namespace AmanziEOS
} // namespace Amanzi
