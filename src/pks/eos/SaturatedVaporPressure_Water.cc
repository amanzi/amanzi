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

#include "SaturatedVaporPressure_Water.hh"

namespace Amanzi {
namespace AmanziEOS {

SaturatedVaporPressure_Water::SaturatedVaporPressure_Water(
    Teuchos::ParameterList& plist)
  : plist_(plist),
    ka0_(16.635764),
    ka_(-6096.9385),
    kb_(-2.7111933e-2),
    kc_(1.673952e-5),
    kd_(2.433502) {};


double SaturatedVaporPressure_Water::Pressure(double T) {
  if (T < 100.0 || T > 373.0) {
    Errors::CutTimeStep msg;
    msg << "Invalid temperature, T = " << T;
    Exceptions::amanzi_throw(msg);
  }
  return 100.0 * exp(ka0_ + ka_/T + (kb_ + kc_*T) * T + kd_*log(T));
}


double SaturatedVaporPressure_Water::DPressureDT(double T) {
  if (T < 100.0 || T > 373.0) {
    std::cout << "Invalid temperature, T = " << T << std::endl;
    Errors::Message msg("Cut time step");
    Exceptions::amanzi_throw(msg);
  }
  return Pressure(T) * (-ka_/(T*T) + kb_ + 2.0*kc_*T + kd_/T);
}

}  // namespace AmanziEOS
}  // namespace Amanzi
