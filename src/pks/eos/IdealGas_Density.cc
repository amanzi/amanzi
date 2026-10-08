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

IdealGas_Density::IdealGas_Density(Teuchos::ParameterList& plist)
  : EOS_Density(plist)
{
  R_ = plist.get<double>("ideal gas constant", 8.31446261815324);
  M_ = plist.get<double>("molar mass");
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

} // namespace AmanziEOS
} // namespace Amanzi
