/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the flow component of the Amanzi code.

*/

#ifndef AMANZI_MULTIPHASE_CONSTANTS_HH_
#define AMANZI_MULTIPHASE_CONSTANTS_HH_

namespace Amanzi {
namespace Multiphase {

// time intervals
const double MULTIPHASE_INITIAL_DT = 1e-8;
const double MULTIPHASE_MAXIMUM_DT = 3.15e+10; // 1000 years

const double MULTIPHASE_WRM_REGULARIZATION_INTERVAL = 0.01;
const double MULTIPHASE_WRM_REGULARIZATION_MAX_GRADIENT = 1.0e+8;
const double MULTIPHASE_WRM_EXCEPTION = -1.0; // triggers exception

const int MULTIPHASE_PHASE_LIQUID = 1; // phases
const int MULTIPHASE_PHASE_GAS = 2;
const int MULTIPHASE_PHASE_NAPL = 3;

} // namespace Multiphase
} // namespace Amanzi

#endif
