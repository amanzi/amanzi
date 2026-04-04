/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Revised Release on the IAPWS-95 formulation
  for the Thermodynamic Properties of Water and Steam.
*/

#ifndef AMANZI_IAPWS95_HH_
#define AMANZI_IAPWS95_HH_

#include <tuple>

#include "Teuchos_ParameterList.hpp"

#include "IAPWS97.hh"

namespace Amanzi {
namespace AmanziEOS {

class IAPWS95 {
 public:
  IAPWS95(Teuchos::ParameterList& plist) : eos97_(plist) {};
  ~IAPWS95() {};

  std::tuple<Properties, Properties, Properties> ThermodynamicsPT(double p, double T);
  std::tuple<Properties, Properties, Properties> ThermodynamicsRhoT(double rho, double T);

  std::array<double, 6> IdealGasPart(double rho, double T);
  std::array<double, 6> ResidualPart(double rho, double T);

  Properties ExtendProperies(double rho, const Properties& prop);

  std::tuple<double, double, double> SaturationLine(double T);
  double VaporPressure(double T);
  double DensityLiquid(double T);
  double DensityVapor(double T);

  // other properties
  double ThermalConductivity(double rho, double T, Properties& prop) {
    return eos97_.ThermalConductivity(rho, T, prop);
  }
  double Viscosity(double rho, double T) { return eos97_.Viscosity(rho, T); }

  // supporting functions
  int get_itrs() { return itrs_; }
  void Print(Properties& prop);

 public:
  int itrs_;

  // static constants
  // IF97 uses different value for R compared to IAPWS95 formulation
  static constexpr double PC = 22.064;     // critical pressure, MPa
  static constexpr double TC = 647.096;    // critical temperature, K
  static constexpr double TT = 273.16;
  static constexpr double RHOC = 322.0;    // critical density, kg/m3
  static constexpr double R = 0.46151805;  // specific gas constant, kJ/kg/K

  // deal-gas part of Helmholtz free energy
  static constexpr double n0[9] = {
    0.0, -8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315, 1.2795, 0.96956, 0.24873
  };
  static constexpr double gamma0[5] = {
     1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105
  };

  // residual part of Helmholtz free energy
  static constexpr int d1[7] = { 1, 1, 1, 2, 2, 3, 4 };
  static constexpr int d2[44] = {
    1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2,  2, 2, 3, 4, 4, 
    4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3,  4, 4, 5, 14, 3, 6, 6, 6
  };
  static constexpr int d3[3] = { 3, 3, 3 };

  static constexpr double t1[7] = { -0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0 };
  static constexpr int t2[44] = {
    4,   6, 12,  1,  5, 4, 2, 13, 9, 3, 4, 11, 4, 13,  1,  7,  1,  9, 10, 10,  3,  7, 
    10, 10,  6, 10, 10, 1, 2,  3, 4, 8, 6,  9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50
  };
  static constexpr int t3[3] = { 0, 1, 4 };

  static constexpr int c2[44] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6
  };

  static constexpr int Nd_max = 16;
  static constexpr int Nt_max = 51;
  static constexpr int Nc_max = 7;

  static constexpr double n1[7] = {
     0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418,
    -0.26145533859358,   -0.78199751687981e-2, 0.88089493102134e-2
  };
  static constexpr double n2[44] = {
    -0.66856572307965,     0.20433810950965,     -0.66212605039687e-4, -0.19232721156002,
    -0.25709043003438,     0.16074868486251,     -0.40092828925807e-1,  0.39343422603254e-6,
    -0.75941377088144e-5,  0.56250979351888e-3,  -0.15608652257135e-4,  0.11537996422951e-8,
     0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9, -0.10793600908932,
     0.17611491008752e-1,  0.22132295167546,     -0.40247669763528,     0.58083399985759,
     0.49969146990806e-2, -0.31358700712549e-1,  -0.74315929710341,     0.47807329915480,
     0.20527940895948e-1, -0.13636435110343,      0.14180634400617e-1,  0.83326504880713e-2,
    -0.29052336009585e-1,  0.38615085574206e-1,  -0.20393486513704e-1, -0.16554050063734e-2,
     0.19955571979541e-2,  0.15870308324157e-3,  -0.16388568342530e-4,  0.43613615723811e-1,
     0.34994005463765e-1, -0.76788197844621e-1,   0.22446277332006e-1, -0.62689710414685e-4,
    -0.55711118565645e-9, -0.19905718354408,      0.31777497330738,    -0.11841182425981
  };
  static constexpr double n3[3] = {
    -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4
  };
  static constexpr double n4[4] = { -0.14874640856724, 0.31806110878444 };

  static constexpr double alpha3[3] = { 20.0, 20.0, 20.0 };
  static constexpr double beta3[3] = { 150.0, 150.0, 250.0 };
  static constexpr double gamma3[3] = { 1.21, 1.21, 1.25 };

  static constexpr double a4[2] = { 3.5, 3.5 };
  static constexpr double b4[2] = { 0.85, 0.95 };
  static constexpr double A[2] = { 0.32, 0.32 };
  static constexpr double B[2] = { 0.2, 0.2 };
  static constexpr double C[2] = { 28.0, 32.0 };
  static constexpr double D[2] = { 700.0, 800.0 };
  static constexpr double beta4[2] = { 0.3, 0.3 };

  // vapor pressure
  static constexpr double PV_n[6] = { 
    -7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502
  };
  static constexpr double PV_k[6] = {
     1.0, 1.5, 3.0, 3.5, 4.0, 7.5
  };

  // density liquid and vapor
  static constexpr double RHOl_n[6] = { 
     1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5
  };
  static constexpr double RHOl_k[6] = { 
     1.0, 2.0, 5.0, 16.0, 43.0, 110.0
  };

  static constexpr double RHOv_n[6] = { 
    -2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581, -63.9201063

  };
  static constexpr double RHOv_k[6] = { 
     2.0, 4.0, 8.0, 18.0, 37.0, 71.0
  };

 private:
  IAPWS97 eos97_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif


