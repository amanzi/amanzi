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
const int FLOW_TIME_INTEGRATION_PICARD = 1;
const int FLOW_TIME_INTEGRATION_BACKWARD_EULER = 2;  // Only for testing.
const int FLOW_TIME_INTEGRATION_BDF1 = 3;
const int FLOW_TIME_INTEGRATION_BDF2 = 4;

// time intervals
const double FLOW_INITIAL_DT = 1e-8;  // [sec]
const double FLOW_MAXIMUM_DT = 3.15e+10;  // [sec] 1000 years
const double FLOW_YEAR = 3.15576e+7;  // [sec]

const double FLOW_PRESSURE_ATMOSPHERIC = 101325.0;

const double FLOW_WRM_VANGENUCHTEN_L = 0.5;
const double FLOW_WRM_BROOKS_COREY_L = 0.5;
const double FLOW_WRM_REGULARIZATION_INTERVAL = 0.0;
const double FLOW_WRM_EXCEPTION = -1.0;  // will trigger exception

const int FLOW_PERMFLAG_NONE = 0;
const int FLOW_PERMFLAG_AVERAGE = 1;
const int FLOW_PERMFLAG_UPWIND = 2;
const int FLOW_PERMFLAG_INTERFACE = 3;
const int FLOW_PERMFLAG_BOUNDARY = 4;

const int FLOW_FV_TPFA = 6;

const int FLOW_TI_ERROR_CONTROL_PRESSURE = 1;  // binary mask for error control
const int FLOW_TI_ERROR_CONTROL_SATURATION = 2;
const int FLOW_TI_ERROR_CONTROL_RESIDUAL = 4;

const double FLOW_TI_ABSOLUTE_TOLERANCE = 1.0;  // defaults for time integrations
const double FLOW_TI_RELATIVE_TOLERANCE = 0.0;
const double FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE = 1e-6;
const int FLOW_TI_MAX_ITERATIONS = 400;

const int FLOW_DT_ADAPTIVE = 1;

const int FLOW_INTERNAL_ERROR = 911;  // contact (lipnikov@lanl.gov)

const int FLOW_UPWIND_UPDATE_TIMESTEP = 1;
const int FLOW_UPWIND_UPDATE_ITERATION = 2;

}  // namespace Flow
}  // namespace Amanzi

#endif

