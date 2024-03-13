/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernel

  Common definitions of for Process Kernels.
*/

#ifndef AMANZI_PK_COMMON_DEFS_HH_
#define AMANZI_PK_COMMON_DEFS_HH_

#include <string>
#include <utility>

namespace Amanzi {
namespace CommonDefs {

const double IDEAL_GAS_CONSTANT_R = 8.314462175;
const double MOLAR_MASS_H2O = 0.0180153333333; // [kg/mol]

// constant properties (T = 293.15 K)
const double ISOTHERMAL_VISCOSITY = 1.002e-3;

const int BOUNDARY_FUNCTION_ACTION_NONE = 0;
const int BOUNDARY_FUNCTION_ACTION_HEAD_RELATIVE = 1;

const int DOMAIN_FUNCTION_ACTION_NONE = 0;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME = 1;
const int DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY = 2;

typedef std::pair<std::string, int> Action;

} // namespace CommonDefs
} // namespace Amanzi

#endif
