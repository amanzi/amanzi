/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  EOS for an ideal gas.
*/

#include "IdealGas_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

IdealGas_Density::IdealGas_Density(Teuchos::ParameterList& eos_plist) : EOS_Density(eos_plist)
{
  InitializeFromPlist_();
};


double
IdealGas_Density::MolarDensity(double T, double p)
{
  return p / (R_ * T);
};


double
IdealGas_Density::DMolarDensityDT(double T, double p)
{
  return -p / (R_ * T * T);
};


double
IdealGas_Density::DMolarDensityDp(double T, double p)
{
  return 1.0 / (R_ * T);
};


void
IdealGas_Density::InitializeFromPlist_()
{
  R_ = 8.31446261815324;
  M_ = eos_plist_.get<double>("molar mass of gas");
};

} // namespace AmanziEOS
} // namespace Amanzi
