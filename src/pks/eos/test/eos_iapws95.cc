/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include "UnitTest++.h"
#include <array>
#include "IAPWS95.hh"

TEST(EOS_IAPWS95)
{
  using namespace Amanzi::AmanziEOS;

  Teuchos::ParameterList plist;
  IAPWS95 eos(plist);

  // ideal gas part
  const auto& g0 = eos.IdealGasPart(838.025, 500.0);
  CHECK_CLOSE(0.204797733e1, g0[0], 1e-8);
  CHECK_CLOSE(0.384236747e0, g0[1], 1e-9);
  CHECK_CLOSE(0.904611106e1, g0[2], 1e-8);
  CHECK_CLOSE(-0.147637878e0, g0[3], 1e-9);
  CHECK_CLOSE(0.0, g0[4], 1e-9);
  CHECK_CLOSE(-0.193249185e1, g0[5], 1e-8);

  // residual part of Helmholtz energy
  auto g = eos.ResidualPart(838.025, 500.0);
  CHECK_CLOSE(-0.342693206e1, g[0], 1e-8);
  CHECK_CLOSE(-0.364366650e0, g[1], 1e-9);
  CHECK_CLOSE(-0.581403435e1, g[2], 1e-8);
  CHECK_CLOSE(0.856063701e0, g[3], 1e-9);
  CHECK_CLOSE(-0.112176915e1, g[4], 1e-8);
  CHECK_CLOSE(-0.223440737e1, g[5], 1e-8);

  g = eos.ResidualPart(358.0, 647.0);
  CHECK_CLOSE(-1.212026565041463, g[0], 1e-12);
  CHECK_CLOSE(-0.714012024371285, g[1], 1e-10);
  CHECK_CLOSE(-3.217225007751656, g[2], 1e-12);
  CHECK_CLOSE(0.475730695645689, g[3], 1e-12);
  CHECK_CLOSE(-1.332147204361430, g[4], 1e-8);
  CHECK_CLOSE(-9.960295065592888, g[5], 1e-8);

  // (rho, T) - variables
  auto [prop, liquid, vapor] = eos.ThermodynamicsRhoT(996.5560, 300.0);
  CHECK_CLOSE(0.0992418350, prop.p, 1e-9);
  CHECK_CLOSE(4.13018112, prop.cv, 1e-8);
  CHECK_CLOSE(1501.51914, prop.w, 1e-5);
  CHECK_CLOSE(0.393062643, prop.s, 1e-9);

  std::tie(prop, liquid, vapor) = eos.ThermodynamicsRhoT(0.435, 500.0);
  CHECK_CLOSE(0.999679423e-1, prop.p, 1e-10);
  CHECK_CLOSE(1.50817541, prop.cv, 1e-8);
  CHECK_CLOSE(548.31425, prop.w, 1e-5);
  CHECK_CLOSE(7.944882714, prop.s, 1e-9);

  std::tie(prop, liquid, vapor) = eos.ThermodynamicsRhoT(358.0, 647.0);
  CHECK_CLOSE(0.220384756e2, prop.p, 1e-7);
  CHECK_CLOSE(0.618315728e1, prop.cv, 1e-8);
  CHECK_CLOSE(0.252145078e3, prop.w, 1e-6);
  CHECK_CLOSE(0.432092307e1, prop.s, 1e-8);

  std::tie(prop, liquid, vapor) = eos.ThermodynamicsRhoT(0.241, 900.0);
  CHECK_CLOSE(0.100062559e0, prop.p, 1e-9);
  CHECK_CLOSE(0.175890657e1, prop.cv, 1e-8);
  CHECK_CLOSE(0.724027147e3, prop.w, 1e-6);
  CHECK_CLOSE(0.916653194e1, prop.s, 1e-8);

  // (p, T) - variables
  std::tie(prop, liquid, vapor) = eos.ThermodynamicsPT(0.992418352e-1, 300.0);
  CHECK_CLOSE(0.9965560e3, prop.rho, 1e-6);
  CHECK_CLOSE(0.413018112e1, prop.cv, 1e-8);
  CHECK_CLOSE(0.150151914e4, prop.w, 1e-5);
  CHECK_CLOSE(0.393062643, prop.s, 1e-9);

  std::tie(prop, liquid, vapor) = eos.ThermodynamicsPT(20.0000690, 900.0);
  CHECK_CLOSE(0.5261500e2, prop.rho, 1e-6);
  CHECK_CLOSE(0.193510526e1, prop.cv, 1e-8);
  CHECK_CLOSE(0.698445674e3, prop.w, 1e-6);
  CHECK_CLOSE(0.659070225e1, prop.s, 1e-8);

  std::tie(prop, liquid, vapor) = eos.ThermodynamicsPT(22.0384756, 647.0);
  CHECK_CLOSE(358.0, prop.rho, 3e-4); // we are close to the critical point

  // two-phase
  double rhol, rhov, ps;
  rhol = eos.DensityLiquid(273.16);
  rhov = eos.DensityVapor(273.16);
  CHECK_CLOSE(999.789, rhol, 1e-3);
  CHECK_CLOSE(0.485426e-2, rhov, 1e-8);

  rhol = eos.DensityLiquid(373.1243);
  rhov = eos.DensityVapor(373.1243);
  CHECK_CLOSE(958.365, rhol, 1e-3);
  CHECK_CLOSE(0.597586, rhov, 1e-6);

  rhol = eos.DensityLiquid(647.096);
  rhov = eos.DensityVapor(647.096);
  CHECK_CLOSE(322.0, rhol, 1e-8);
  CHECK_CLOSE(322.0, rhov, 1e-8);

  std::tie(rhol, rhov, ps) = eos.SaturationLine(275.0);
  CHECK_CLOSE(0.999887406e3, rhol, 1e-6);
  CHECK_CLOSE(0.550664919e-2, rhov, 1e-11);
  CHECK_CLOSE(0.698451167e-3, ps, 1e-12);

  std::tie(rhol, rhov, ps) = eos.SaturationLine(450.0);
  CHECK_CLOSE(0.890341250e3, rhol, 1e-6);
  CHECK_CLOSE(0.481200360e1, rhov, 1e-8);
  CHECK_CLOSE(0.932203564, ps, 1e-9);

  std::tie(rhol, rhov, ps) = eos.SaturationLine(625.0);
  CHECK_CLOSE(0.567090385e3, rhol, 1e-6);
  CHECK_CLOSE(0.118290280e3, rhov, 1e-6);
  CHECK_CLOSE(0.169082693e2, ps, 1e-7);
}
