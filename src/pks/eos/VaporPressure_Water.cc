/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Saturated vapor pressure for vapor over water or ice, Sonntag (1990)
  Note that in his formula, the saturated vapor pressure is defined in hPa.
*/

#include <cmath>
#include "errors.hh"
#include "VaporPressure_Water.hh"

namespace Amanzi {
namespace AmanziEOS {

VaporPressure_Water::VaporPressure_Water(Teuchos::ParameterList& plist) :
    plist_(plist),
    ka0_(16.635764),
    ka_(-6096.9385),
    kb_(-2.7111933e-2),
    kc_(1.673952e-5),
    kd_(2.433502) {}


double VaporPressure_Water::SaturatedVaporPressure(double T) {
  if (T < 100.0 || T > 373.0) {
    std::cout << "Invalid temperature, T = " << T << std::endl;
    Exceptions::amanzi_throw(Errors::CutTimeStep());
  }
  return 100.0 * exp(ka0_ + ka_/T + (kb_ + kc_*T) * T + kd_*log(T));
}


double VaporPressure_Water::DSaturatedVaporPressureDT(double T) {
  if (T < 100.0 || T > 373.0) {
    std::cout << "Invalid temperature, T = " << T << std::endl;
    Errors::Message msg("Cut time step");
    Exceptions::amanzi_throw(msg);
  }
  return SaturatedVaporPressure(T) * (-ka_/(T*T) + kb_ + 2.0*kc_*T + kd_/T);
}

}  // namespace AmanziEOS
}  // namespace Amanzi
