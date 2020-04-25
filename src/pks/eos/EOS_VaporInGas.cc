/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for a combination of air and vapor pressure. Mass density is not
  available, not because it can't be calculated, but because it depends 
  upon omega. It's not really needed, and if it were, would not fit the 
*/

#include "EOSFactory.hh"
#include "EOS_VaporInGas.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor takes a parameter list with EOS parameters.
******************************************************************* */
EOS_VaporInGas::EOS_VaporInGas(Teuchos::ParameterList& eos_plist) :
    eos_plist_(eos_plist) {
  InitializeFromPlist_();
}


double EOS_VaporInGas::MolarDensity(double T, double p) {
  return gas_eos_->MolarDensity(T,p);
}


double EOS_VaporInGas::DMolarDensityDT(double T, double p) {
  return gas_eos_->DMolarDensityDT(T,p);
}


double EOS_VaporInGas::DMolarDensityDp(double T, double p) {
  return gas_eos_->DMolarDensityDp(T,p);
}


void EOS_VaporInGas::InitializeFromPlist_() {
  Teuchos::ParameterList gas_eos_plist = eos_plist_.sublist("gas EOS parameters");
  EOSFactory eos_factory;
  gas_eos_ = eos_factory.CreateEOS(gas_eos_plist);
}

}  // namespace AmanziEOS
}  // namespace Amanzi
