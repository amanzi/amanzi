/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_IAPWS_BASE_HH_
#define AMANZI_IAPWS_BASE_HH_

namespace Amanzi {
namespace AmanziEOS {

enum class Phase_t : int {
  None = 0,
  CompressibleLiquid = 1,
  Gas = 2,
  CriticalPoint = 3,
  SaturatedVapor = 4,
  SaturatedLiquid = 5,
  TwoPhases = 6,
  SupercriticalLiquid = 7,
  Vapor = 8,
  Liquid = 9
};


struct Properties {
  double p = 0.0; // pressure, MPa
  double T = 0.0; // temeparture, K
  double rho = 0.0;
  double v = 0.0; // specific volume, m3/kg
  double h = 0.0; // specific enthalpy, kJ/kg
  double u = 0.0; // specific internal energy, kJ/kg
  double s = 0.0; // specific entropy, kJ/kg/K
  double cp = 0.0; // specific isobaric heat capacity, kJ/kg/K
  double cv = 0.0; // specific isocoric heat capacity, kJ/kg/K
  double w = 0.0; // speed of sound, m/s
  double kt = 0.0; // isothermal compressibility, 1/MPa
  double av = 0.0; // isobaric cubic expansion coefficient, 1/K
  double ap = 0.0; // relative pressure coefficient, 1/K
  double bp = 0.0; // isothermal stress coefficient, kg/m3

  double helmholtz = 0.0; // specific Helmholtz free energy, kJ/kg
  double gibbs = 0.0; // specific Gibbs free energy, kJ/kg

  double mu = 0.0; // dynamic viscosity, Pa s
  double k = 0.0; // thermal conductivity, W/m/K

  double sigma = 0.0; // surface tension, N/m

  Phase_t phase = Phase_t::None; // phase id
  double x = 0.0; // vapor quality

  int rgn = 0; // region id
};


class IAPWS {
 public:
  IAPWS(Teuchos::ParameterList& plist) {
    k_enhancement_ = plist.get<bool>("thermal conductivity enhancement", false);
    mu_enhancement_ = plist.get<bool>("viscosity enhancement", false);
  }
  ~IAPWS() {};

  virtual std::tuple<Properties, Properties, Properties> ThermodynamicsPT(double p, double T) = 0;
  virtual std::tuple<Properties, Properties, Properties> ThermodynamicsPH(double p, double h) = 0;

  double ThermalConductivity(double rho, double T, Properties& prop);
  double ThermalConductivityBase(double rho, double T) const;
  double ThermalConductivityCriticalEnhancement(const Properties& prop) const;

  double SurfaceTension(double T);

  double Viscosity(double rho, double T, const Properties& prop);
  double ViscosityBase(double rho, double T) const;
  double ViscosityCriticalEnhancement(const Properties& prop) const;
  double ViscosityBaseDerivativeRho(double rho, double T) const;
  double ViscosityBaseDerivativeT(double rho, double T) const;

 public:
  static constexpr double TC = 647.096;  // critical temperature, K
  static constexpr double PC = 22.064;   // critical pressure, MPa
  static constexpr double RHOC = 322.0;  // critical density, kg/m3

  // thermal conductivity
  static constexpr int N0_ThCond = 5;
  static constexpr double thermal_cond_n0[N0_ThCond] = {
    2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4
  };

  static constexpr int N1_ThCond = 6;
  static constexpr double thermal_cond_n1[N0_ThCond][N1_ThCond] = {
    { 1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634,  0.00609859258 },
    { 2.33771842, -2.78843778,  1.53616167, -0.463045512,  0.0832827019, -0.00719201245 },
    { 2.19650529, -4.54580785,  3.55777244, -1.40944978,   0.275418278,  -0.0205938816 },
    {-1.21051378,  1.60812989, -0.621178141, 0.0716373224, 0.0,           0.0 },
    {-2.7203370,   4.57586331, -3.18369245,  1.1168348,   -0.19268305,    0.012913842 }
  };

  // -- critical enhancement constants
  static constexpr double Lambda = 177.8514;
  static constexpr double qD_inv_k = 0.40e-9;  // m
  static constexpr double xi0_k = 0.13e-9;  // m

  static constexpr double nu = 0.630;
  static constexpr double gamma  = 1.239;
  static constexpr double Gamma0 = 0.06;

  static constexpr int M2_ThCond = 6;
  static constexpr int N2_ThCond = 5;
  static constexpr double thermal_cond_zeta[M2_ThCond][N2_ThCond] = {
    {  6.53786807199516,  6.52717759281799,  5.35500529896124,  1.55225959906681,  1.11999926419994 },
    { -5.61149954923348, -6.30816983387575, -3.96415689925446,  0.464621290821181, 0.595748562571649 },
    {  3.39624167361325,  8.08379285492595,  8.91990208918795,  8.93237374861479,  9.88952565078920 },
    { -2.27492629730878, -9.82240510197603,-12.0338729505790, -11.0321960061126, -10.3255051147040 },
    { 10.2631854662709,  12.1358413791395,   9.19494865194302,  6.16780999933360,  4.66861294457414 },
    {  1.97815050331519, -5.54349664571295, -2.16866274479712, -0.965458722086812,-0.503243546373828 },
  };
  static constexpr double thermal_cond_rho[N2_ThCond] = {
    0.310559006, 0.776397516, 1.242236025, 1.863354037, 100.0
  };

  // dynamic viscosity 
  static constexpr int N0_Visc = 4;
  static constexpr double viscosity_n0[N0_Visc] = {
    1.67752, 2.20462, 0.6366564, -0.241605
  };

  static constexpr int N1_Visc = 21;
  static constexpr int viscosity_k1[N1_Visc] = {
    0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5
  };
  static constexpr int N1_Visc_kmax = 6;
  static constexpr int viscosity_l1[N1_Visc] = {
    0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6
  };
  static constexpr int N1_Visc_lmax = 7;
  static constexpr double viscosity_n1[N1_Visc] = {
    0.520094,     0.850895e-1, -0.108374e1, -0.289555,  0.222531,    0.999115,
    0.188797e1,   0.126613e1,   0.120573,   -0.281378, -0.906851,   -0.772479,
   -0.489837,    -0.257040,     0.161913,    0.257399, -0.325372e-1, 0.698452e-1,
    0.872102e-2, -0.435673e-2, -0.593264e-3
  };

  // other data
  static constexpr double x_mu = 0.068;

  static constexpr double qC_inv_mu = 1.9;  // reciprocal wave numbers, nm
  static constexpr double qD_inv_mu = 1.1;  // nm

  static constexpr double xi0_mu = 0.13;  // nm
  static constexpr double xi_switch = 0.3817016416;  // nm

 private:
  double ThermalConductivityZetaFunction_(double rhor) const;

 protected:
  bool k_enhancement_;
  bool mu_enhancement_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
