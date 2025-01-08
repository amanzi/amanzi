/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  EOS

  Viscosity for ideal gas based on the Sutherland law.
*/

#include "errors.hh"
#include "IdealGas_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

IdealGas_Viscosity::IdealGas_Viscosity(Teuchos::ParameterList& plist)
  : EOS_Viscosity(plist), T0_(273.0)
{
  mu0_ = plist.get<double>("reference viscosity");
  T0_ = plist.get<double>("reference temperature");
  S0_ = plist.get<double>("Sutherland constant");

  factor_ = mu0_ * (T0_ + S0_) / std::pow(T0_, 1.5);
}


double
IdealGas_Viscosity::Viscosity(double T, double p)
{
  return factor_ * std::pow(T, 1.5) / (T + S0_);
}


double
IdealGas_Viscosity::DViscosityDT(double T, double p)
{
  double tmp = T + S0_;
  return factor_ * std::sqrt(T) * (1.5 * tmp - T) / tmp / tmp;
}

} // namespace AmanziEOS
} // namespace Amanzi
