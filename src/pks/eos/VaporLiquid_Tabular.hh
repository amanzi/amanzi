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

static std::map<std::string, std::vector<double>> TempRange = { { "CO2", { 274.19, 642.66 } },
                                                                { "Xe", { 273.22, 543.36 } },
                                                                { "CH4", { 275.44, 473.46 } } };

static std::map<std::string, std::vector<double>> EquilibriumCoef = {
  { "CO2", { 1672.9376, 28.1751, -112.4619, 85.3807 } },
  { "Xe", { 2022, 8375, 16.7913, -61.2401, 41.9236 } },
  { "CH4", { 2215, 6977, -0.1089, -6.6240, 4.6789 } }
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
