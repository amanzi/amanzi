/*
  Shallow Water PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
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
  hydrostatic_pressure_force_type_ = plist.get<std::string>("hydrostatic pressure force type", "shallow water");
  pipe_diameter_ = plist.get<double>("pipe diameter", 1.0);
  celerity_ = plist.get<double>("celerity", 100);
}


/* ******************************************************************
* Numerical flux
****************************************************************** */
std::vector<double>
NumericalFlux_Rusanov::Compute(const std::vector<double>& UL, const std::vector<double>& UR)
{
  std::vector<double> FL, FR, F(3);

  double hL, uL, vL, hR, uR, vR, factor;
  double eps = 1.e-6;

  // SW conservative variables: (h, hu, hv)

  hL = UL[0];
  factor = 2.0 * hL / (hL * hL + std::fmax(hL * hL, eps * eps));
  uL = factor * UL[1];
  vL = factor * UL[2];

  hR = UR[0];
  factor = 2.0 * hR / (hR * hR + std::fmax(hR * hR, eps * eps));
  uR = factor * UR[1];
  vR = factor * UR[2];

  FL = PhysicalFlux(UL);
  FR = PhysicalFlux(UR);

  double SL, SR, Smax;

  SL = std::fabs(uL) + std::sqrt(g_ * hL);
  SR = std::fabs(uR) + std::sqrt(g_ * hR);

  Smax = std::max(SL, SR);
  lambda_max_ = Smax; // FIXME, we probably need only max speed
  lambda_min_ = -Smax;

  for (int i = 0; i < 3; i++) { F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * Smax * (UR[i] - UL[i]); }

  return F;
}

} // namespace ShallowWater
} // namespace Amanzi
