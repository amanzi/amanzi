/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Tabulated vapor-Liquid distribution based on J. Phys. Chem. Ref. Data, 
  Vol.32, number 2, 2003
   
   kD(T) = xg / xl 

  where xg and xl are gas and liquid more fractions in equilibrium.
*/

#ifndef AMANZI_EOS_VAPOR_LIQUID_TABULAR_HH_
#define AMANZI_EOS_VAPOR_LIQUID_TABULAR_HH_

#include <map>

#include "dbc.hh"
#include "Factory.hh"

#include "VaporLiquid.hh"

namespace Amanzi {
namespace AmanziEOS {

static double kB1 = 1.99274064;
static double kB2 = 1.09965342;
static double kB3 = -0.510839303;
static double kB4 = -1.75493479;
static double kB5 = -45.5170352;
static double kB6 = -6.74694450e+5;

static double kq = -0.023767;
static double kTc = 647.096;

static std::map<std::string, std::vector<double>> TempRange = {
  { "He", { 273.21, 553.18 } },   { "Ne", { 273.20, 543.36 } },  { "Ar", { 273.19, 568.36 } },
  { "Kr", { 273.19, 525.56 } },   { "Xe", { 273.22, 574.85 } },  { "H2", { 273.15, 636.09 } },
  { "N2", { 278.12, 636.46 } },   { "O2", { 274.15, 616.52 } },  { "CO", { 278.15, 588.67 } },
  { "CO2", { 274.19, 642.66 } },  { "H2S", { 273.15, 533.09 } }, { "CH4", { 275.46, 633.11 } },
  { "C2H6", { 275.44, 473.46 } }, { "SF6", { 283.14, 505.55 } }
};

static std::map<std::string, std::vector<double>> EquilibriumCoef = {
  { "He", { 2267.4082, -2.9616, -3.2604, 7.8819 } },
  { "Ne", { 2507.3022, -38.6955, 110.3992, -71.9096 } },
  { "Ar", { 2310.5463, -46.7034, 160.4066, -118.3043 } },
  { "Kr", { 2276.9722, -61.1494, 214.0117, -159.0407 } },
  { "Xe", { 2022.8375, 16.7913, -61.2401, 41.9236 } },
  { "H2", { 2286.4159, 11.3397, -70.7279, 63.0631 } },
  { "N2", { 2388.8777, -14.9593, 42.0179, -29.4396 } },
  { "O2", { 2305.0674, -11.3240, 25.3224, -15.6449 } },
  { "CO", { 2346.2291, -57.6317, 204.5324, -152.6377 } },
  { "CO2", { 1672.9376, 28.1751, -112.4619, 85.3807 } },
  { "H2S", { 1319.1205, 14.1571, -46.8361, 33.2266 } },
  { "CH4", { 2215.6977, -0.1089, -6.6240, 4.6789 } },
  { "C2H6", { 2143.8121, 6.8859, -12.6084, 0.0000 } },
  { "SF6", { 2871.7265, -66.7556, 229.7191, -172.7400 } }
};


class VaporLiquid_Tabular : public VaporLiquid {
 public:
  VaporLiquid_Tabular(const std::string& name);

  virtual double k(double T) const;

  virtual double DkDT(double T) const { return 0.0; } // FIXME

 private:
  double E_, F_, G_, H_, Tmin_, Tmax_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
