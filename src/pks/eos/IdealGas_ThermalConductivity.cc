/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Thermal conductivity of ideal gas based on the Sutherland law.
*/

#include "errors.hh"
#include "IdealGas_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

IdealGas_ThermalConductivity::IdealGas_ThermalConductivity(Teuchos::ParameterList& plist)
  : EOS_ThermalConductivity(plist), T0_(273.0)
{
  k0_ = plist.get<double>("reference conductivity");
  T0_ = plist.get<double>("reference temperature");
  S0_ = plist.get<double>("Sutherland constant");

  factor_ = k0_ * (T0_ + S0_) / std::pow(T0_, 1.5);
}


double
IdealGas_ThermalConductivity::ThermalConductivity(double p, double T, double pgi)
{
  return factor_ * std::pow(T, 1.5) / (T + S0_);
}


double
IdealGas_ThermalConductivity::DThermalConductivityDT(double p, double T, double phi)
{
  double tmp = T + S0_;
  return factor_ * std::sqrt(T) * (1.5 * tmp - T) / tmp / tmp;
}

} // namespace AmanziEOS
} // namespace Amanzi
