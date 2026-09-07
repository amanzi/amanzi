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
  virtual std::tuple<Properties, Properties, Properties> ThermodynamicsPT(double p, double T) = 0;
  virtual std::tuple<Properties, Properties, Properties> ThermodynamicsPH(double p, double h) = 0;

  virtual double ThermalConductivity(double rho, double T, Properties& prop) = 0;
  virtual double Viscosity(double rho, double T) = 0;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
