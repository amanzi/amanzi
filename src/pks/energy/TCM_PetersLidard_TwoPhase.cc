/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon 

  Linear interpolant of thermal conductivity.
*/

#include <cmath>
#include "TCM_PetersLidard_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

TCM_PetersLidard_TwoPhase::TCM_PetersLidard_TwoPhase(Teuchos::ParameterList& plist) : plist_(plist)
{
  InitializeFromPlist_();
}


double
TCM_PetersLidard_TwoPhase::ThermalConductivity(double poro, double sat_liq)
{
  double k_dry = (d_ * (1.0 - poro) * k_rock_ + k_gas_ * poro) / (d_ * (1.0 - poro) + poro);
  double k_sat = pow(k_rock_, (1.0 - poro)) * pow(k_liquid_, poro);
  double kersten = pow(sat_liq + eps_, alpha_);
  return k_dry + (k_sat - k_dry) * kersten;
}


void
TCM_PetersLidard_TwoPhase::InitializeFromPlist_()
{
  d_ = 0.053; // unitless empericial parameter

  eps_ = plist_.get<double>("epsilon", 1.0e-10);
  alpha_ = plist_.get<double>("unsaturated alpha");
  k_rock_ = plist_.get<double>("thermal conductivity of rock");
  k_liquid_ = plist_.get<double>("thermal conductivity of liquid", 0.6065);
  k_gas_ = plist_.get<double>("thermal conductivity of gas");
}

} // namespace Energy
} // namespace Amanzi
