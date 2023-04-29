/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS
*/

#include "VaporLiquid_Tabular.hh"

namespace Amanzi {
namespace AmanziEOS {

VaporLiquid_Tabular::VaporLiquid_Tabular(const std::string& name)
{
  E_ = EquilibriumCoef[name][0];
  F_ = EquilibriumCoef[name][1];
  G_ = EquilibriumCoef[name][2];
  H_ = EquilibriumCoef[name][3];

  Tmin_ = TempRange[name][0];
  Tmax_ = TempRange[name][1];
}


double VaporLiquid_Tabular::k(double T) const
{
  double s, tau, f, kD;
  s = T / kTc;
  tau = 1 - s;
  f = kB1 * std::pow(tau, 1.0 / 3)
    + kB2 * std::pow(tau, 2.0 / 3)
    + kB3 * std::pow(tau, 5.0 / 3)
    + kB4 * std::pow(tau, 16.0 / 3)
    + kB5 * std::pow(tau, 43.0 / 3)
    + kB6 * std::pow(tau, 110.0 / 3);

  kD = kq * F_ + E_ * f / T + F_ + G_ * std::pow(tau, 2.0 / 3) + H_ * tau;
  return kD;
}

} // namespace AmanziEOS
} // namespace Amanzi

