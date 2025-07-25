/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Shallow Water PK

*/

#include <cmath>
#include <string>

#include "NumericalFlux_Rusanov.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Numerical flux
****************************************************************** */
NumericalFlux_Rusanov::NumericalFlux_Rusanov(Teuchos::ParameterList& plist)
{
  g_ = plist.get<double>("gravity");
}


/* ******************************************************************
* Numerical flux
****************************************************************** */
std::vector<double>
NumericalFlux_Rusanov::Compute(const std::vector<double>& UL,
                               const std::vector<double>& UR,
                               const double& HPFL,
                               const double& HPFR)
{
  std::vector<double> FL, FR, F(3);

  double hL, uL, hR, uR, factor;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  hL = UL[0];
  factor = 2.0 * hL / (hL * hL + std::fmax(hL * hL, eps * eps));
  uL = factor * UL[1];

  hR = UR[0];
  factor = 2.0 * hR / (hR * hR + std::fmax(hR * hR, eps * eps));
  uR = factor * UR[1];

  FL = PhysicalFlux(UL, HPFL);
  FR = PhysicalFlux(UR, HPFR);

  double SL, SR, Smax;

  SL = std::fabs(uL) + std::sqrt(g_ * hL);
  SR = std::fabs(uR) + std::sqrt(g_ * hR);

  Smax = std::max(SL, SR);
  lambda_max_ = Smax; // FIXME, we probably need only max speed
  lambda_min_ = -Smax;

  for (int i = 0; i < 3; i++) {
    F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * Smax * (UR[i] - UL[i]);
  }

  return F;
}

} // namespace ShallowWater
} // namespace Amanzi
