/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_TIME_INTEGRATION_CONSTANTS_HH_
#define AMANZI_TIME_INTEGRATION_CONSTANTS_HH_


namespace Amanzi {

const double DT_CONTROLLER_ADAPTIVE_INCREASE = 4.0;
const double DT_CONTROLLER_ADAPTIVE_REDUCTION = 0.1;
const double DT_CONTROLLER_ADAPTIVE_SAFETY_FACTOR = 0.9;
const double DT_CONTROLLER_ADAPTIVE_ERROR_TOLERANCE = 1e-10;

} // namespace Amanzi

#endif
