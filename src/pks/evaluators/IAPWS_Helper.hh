/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators
*/

#ifndef AMANZI_EVALUATORS_IAPWS_HELPER_HH_
#define AMANZI_EVALUATORS_IAPWS_HELPER_HH_

#include <cmath>
#include <string>
#include <vector>

namespace Amanzi {
namespace Evaluators {

int constexpr TSPT_t_size = 17;
const std::vector<std::string> TSPT_names = {
  "mass_density", "enthalpy", "specific_volume", "isobaric_heat_capacity",
  "isocoric_heat_capacity", "isothermal_compressibility", "isobaric_expansion_coef", 
  "relative_pressure_coef", "isothermal_stress_coef", "thermal_conductivity", "viscosity",
  "drhodp", "drhodT", "dudp", "dudT", "vv", "vapor_quality"
};

enum class TSPT_t : int {
  RHO = 0,
  H = 1,
  V = 2,
  CP = 3,
  CV = 4,
  KT = 5,
  AV = 6,
  AP = 7,
  BP = 8,
  K = 9,
  MU = 10,
  dRHOdP = 11, // derivatives
  dRHOdT = 12,
  dUdP = 13,
  dUdT = 14,
  VV = 15, // two-phase
  X = 16
};


int constexpr TSPH_t_size = 18;
const std::vector<std::string> TSPH_names = {
  "region", "temperature", "mass_density", "specific_volume", "isobaric_heat_capacity",
  "isocoric_heat_capacity", "isothermal_compressibility", "isobaric_expansion_coef", 
  "relative_pressure_coef", "isothermal_stress_coef", "thermal_conductivity", "viscosity",
  "drhodp", "drhodh", "dtdp", "dtdh",
  "vv", "vapor_quality"
};

enum class TSPH_t : int {
  RGN = 0,
  T = 1,
  RHO = 2,
  V = 3,
  CP = 4,
  CV = 5,
  KT = 6,
  AV = 7,
  AP = 8,
  BP = 9,
  K = 10,
  MU = 11,
  dRHOdP = 12, // derivatives
  dRHOdH = 13,
  dTdP = 14,
  dTdH = 15,
  VV = 16, // two-phase
  X = 17
};

} // namespace Evaluators
} // namespace Amanzi

#endif
