/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for an ideal gas.
*/

#include "EOS_DensityIdealGas.hh"

namespace Amanzi {
namespace AmanziEOS {

EOS_DensityIdealGas::EOS_DensityIdealGas(Teuchos::ParameterList& eos_plist)
  : EOS_Density(eos_plist) {
  InitializeFromPlist_();
};


double EOS_DensityIdealGas::MolarDensity(double T, double p) {
  return p / (R_ * T);
};


double EOS_DensityIdealGas::DMolarDensityDT(double T, double p) {
  return -p / (R_ * T * T);
};


double EOS_DensityIdealGas::DMolarDensityDp(double T, double p) {
  return 1.0 / (R_ * T);
};


void EOS_DensityIdealGas::InitializeFromPlist_()
{
  R_ = 8.31446261815324;
  M_ = eos_plist_.get<double>("molar mass of gas");
};
 
}  // namespace AmanziEOS
}  // namespace Amanzi
