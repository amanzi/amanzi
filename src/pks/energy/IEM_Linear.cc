/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Linear internal energy model -- function of cv and temperature

   IE = cv * (T - Tref)

  UNITS: J/{mol/kg}
*/

#include "IEM_Linear.hh"

namespace Amanzi {
namespace Energy {

IEM_Linear::IEM_Linear(Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
}


double IEM_Linear::InternalEnergy(double temp) {
  return cv_ * (temp - Tref_);
}


void IEM_Linear::InitializeFromPlist_()
{
  if (plist_.isParameter("heat capacity [J/kg-K]")) {
    cv_ = plist_.get<double>("heat capacity [J/kg-K]");
    molar_basis_ = false;
  } else {
    cv_ = plist_.get<double>("heat capacity [J/mol-K]");
    molar_basis_ = true;
  }

  Tref_ = plist_.get<double>("reference temperature", 273.15);
}

}  // namespace Energy
}  // namespace Amanzi
