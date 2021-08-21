/*
  Shallow Water PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
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
std::vector<double> NumericalFlux_CentralUpwind::Compute(
    const std::vector<double>& UL, const std::vector<double>& UR)
{
  std::vector<double> FL, FR, F(3), U_star(3), dU(3);

  double hL, uL, vL, hR, uR, vR, qxL, qyL, qxR, qyR;
  double apx, amx;
  double eps = 1.0e-6, eps_1 = 1.0e-14;

  // SW conservative variables: (h, hu, hv)

  hL  = UL[0];
  qxL = UL[1];
  qyL = UL[2];
  uL  = 2.*hL*qxL/(hL*hL + std::fmax(hL*hL, eps*eps));
  vL  = 2.*hL*qyL/(hL*hL + std::fmax(hL*hL, eps*eps));

  hR  = UR[0];
  qxR = UR[1];
  qyR = UR[2];
  uR  = 2.*hR*qxR/(hR*hR + std::fmax(hR*hR, eps*eps));
  vR  = 2.*hR*qyR/(hR*hR + std::fmax(hR*hR, eps*eps));

  apx = std::max(std::max(uL + std::sqrt(g_*hL), uR + std::sqrt(g_*hR)), 0.0);
  lambda_max_ = apx;

  amx = std::min(std::min(uL - std::sqrt(g_*hL), uR-std::sqrt(g_*hR)), 0.0);
  lambda_min_ = amx;

  FL = PhysicalFlux(UL);
  FR = PhysicalFlux(UR);
  
  for (int i = 0; i < 3; i++) {
    U_star[i] = (apx*UR[i] - amx*UL[i] - (FR[i] - FL[i])) / (apx - amx + eps_1);
  }

  for (int i = 0; i < 3; i++) {
    dU[i] = minmod(UR[i]-U_star[i],U_star[i]-UL[i]);
  }

  for (int i = 0; i < 3; i++) {
    F[i] = (apx*FL[i] - amx*FR[i] + apx*amx*(UR[i] - UL[i] - dU[i])) / (apx - amx + eps_1);
  }

  return F;
}

}  // namespace ShallowWater
}  // namespace Amanzi

