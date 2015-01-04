/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  EOS for an ideal gas (does not implement viscosity at this point!)
*/

#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace Relations {

EOSVaporInGas::EOSVaporInGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {
  InitializeFromPlist_();
}


double EOSVaporInGas::MolarDensity(double T, double p) {
  return gas_eos_->MolarDensity(T,p);
};


double EOSVaporInGas::DMolarDensityDT(double T, double p) {
  return gas_eos_->DMolarDensityDT(T,p);
};


double EOSVaporInGas::DMolarDensityDp(double T, double p) {
  return gas_eos_->DMolarDensityDp(T,p);
};


void EOSVaporInGas::InitializeFromPlist_() {
  Teuchos::ParameterList gas_eos_plist = eos_plist_.sublist("gas EOS parameters");
  EOSFactory eos_factory;
  gas_eos_ = eos_factory.createEOS(gas_eos_plist);
};

}  // namespace Relations
}  // namespace Amanzi
