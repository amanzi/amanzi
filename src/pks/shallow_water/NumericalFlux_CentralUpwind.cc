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

#include "NumericalFlux_CentralUpwind.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Numerical flux
****************************************************************** */
NumericalFlux_CentralUpwind::NumericalFlux_CentralUpwind(Teuchos::ParameterList& plist)
{
  g_ = plist.get<double>("gravity");
}


/* ******************************************************************
* Numerical flux
****************************************************************** */
std::vector<double>
NumericalFlux_CentralUpwind::Compute(const std::vector<double>& UL,
                                     const std::vector<double>& UR,
                                     const double& HPFL,
                                     const double& HPFR)
{
  std::vector<double> FL, FR, F(3), U_star(3), dU(3);

  double hL, uL, hR, uR, qxL, qxR;
  double apx, amx, factor, ghL, ghR;
  double eps0 = 1.0e-12, eps1 = 1.0e-14;

  // SW conservative variables: (h, hu, hv)

  hL = UL[0];
  qxL = UL[1];
  factor = 2 * hL / (hL * hL + std::fmax(hL * hL, eps0));
  uL = qxL * factor;

  hR = UR[0];
  qxR = UR[1];
  factor = 2 * hR / (hR * hR + std::fmax(hR * hR, eps0));
  uR = qxR * factor;

  ghL = std::sqrt(g_ * hL);
  ghR = std::sqrt(g_ * hR);
  apx = std::max(std::max(uL + ghL, uR + ghR), 0.0);
  lambda_max_ = apx;

  amx = std::min(std::min(uL - ghL, uR - ghR), 0.0);
  lambda_min_ = amx;

  FL = PhysicalFlux(UL, HPFL);
  FR = PhysicalFlux(UR, HPFR);

  for (int i = 0; i < 3; i++) {
    U_star[i] = (apx * UR[i] - amx * UL[i] - (FR[i] - FL[i])) / (apx - amx + eps1);
  }

  for (int i = 0; i < 3; i++) {
    dU[i] = minmod(UR[i] - U_star[i], U_star[i] - UL[i]);
  }

  for (int i = 0; i < 3; i++) {
    F[i] = (apx * FL[i] - amx * FR[i] + apx * amx * (UR[i] - UL[i] - dU[i])) / (apx - amx + eps1);
  }

  return F;
}

} // namespace ShallowWater
} // namespace Amanzi
