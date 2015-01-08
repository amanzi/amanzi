/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Linear internal energy model -- function of Cv and temperature

   IE = Cv * (T - T_ref)

  UNITS: J/{mol/kg}
*/

#include "iem_linear.hh"

namespace Amanzi {
namespace Energy {

IEMLinear::IEMLinear(Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
};


double IEMLinear::InternalEnergy(double temp) {
  return Cv_ * (temp - T_ref_);
};


void IEMLinear::InitializeFromPlist_()
{
  if (plist_.isParameter("heat capacity [J/kg-K]")) {
    Cv_ = plist_.get<double>("heat capacity [J/kg-K]");
    molar_basis_ = false;
  } else {
    Cv_ = plist_.get<double>("heat capacity [J/mol-K]");
    molar_basis_ = true;
  }

  T_ref_ = plist_.get<double>("Reference temperature [K]", 273.15);
};

}  // namespace Energy
}  // namespace Amanzi
