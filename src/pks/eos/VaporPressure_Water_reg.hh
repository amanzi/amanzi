/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Saturated vapor pressure for vapor over water or ice, Sonntag (1990)
*/

#include <cmath>
#include "errors.hh"
#include "VaporPressure_Water.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<VaporPressure_Base, VaporPressure_Water> VaporPressure_Water::factory_("water vapor over water/ice");

}  // namespace AmanziEOS
}  // namespace Amanzi
