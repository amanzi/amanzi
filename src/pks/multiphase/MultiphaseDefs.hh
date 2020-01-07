/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTIPHASE_CONSTANTS_HH_
#define AMANZI_MULTIPHASE_CONSTANTS_HH_

namespace Amanzi {
namespace Multiphase {

// special bits for submodels
const int MULTIPHASE_TIME_INTEGRATION_PICARD = 1;
const int MULTIPHASE_TIME_INTEGRATION_BACKWARD_EULER = 2;  // Only for testing.
const int MULTIPHASE_TIME_INTEGRATION_BDF1 = 3;

// time intervals
const double MULTIPHASE_INITIAL_DT = 1e-8;  // [sec]
const double MULTIPHASE_MAXIMUM_DT = 3.15e+10;  // [sec] 1000 years

const double MULTIPHASE_WRM_REGULARIZATION_INTERVAL = 0.0;
const double MULTIPHASE_WRM_EXCEPTION = -1.0;  // will trigger exception

const int MULTIPHASE_TI_ERROR_CONTROL_PRESSURE = 1;  // binary mask for error control
const int MULTIPHASE_TI_ERROR_CONTROL_SATURATION = 2;
const int MULTIPHASE_TI_ERROR_CONTROL_RESIDUAL = 4;

const double MULTIPHASE_TI_NONLINEAR_RESIDUAL_TOLERANCE = 1e-6;  // defaults for time integrations
const int MULTIPHASE_TI_MAX_ITERATIONS = 400;

const int MULTIPHASE_DT_ADAPTIVE = 1;

}  // namespace Flow
}  // namespace Amanzi

#endif

