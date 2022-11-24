/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_CONSTANTS_HH_
#define AMANZI_FLOW_CONSTANTS_HH_

namespace Amanzi {
namespace Flow {

// default parameters for some seepage face models
const double FLOW_BC_SEEPAGE_FACE_IMPEDANCE = 1e-10;       // [sec / m]
const double FLOW_BC_SEEPAGE_FACE_REGULARIZATION = 5000.0; // [Pa]

// time intervals
const double FLOW_INITIAL_DT = 1e-8;     // [sec]
const double FLOW_MAXIMUM_DT = 3.15e+10; // [sec] 1000 years
const double FLOW_YEAR = 3.15576e+7;     // [sec]

const double FLOW_PRESSURE_ATMOSPHERIC = 101325.0;

const double FLOW_WRM_VANGENUCHTEN_L = 0.5;
const double FLOW_WRM_BROOKS_COREY_L = 0.5;
const double FLOW_WRM_REGULARIZATION_INTERVAL = 0.0;
const double FLOW_WRM_EXCEPTION = -1.0; // will trigger exception

// time integration
const int FLOW_TI_ERROR_CONTROL_PRESSURE = 1; // binary mask for error control
const int FLOW_TI_ERROR_CONTROL_SATURATION = 2;
const int FLOW_TI_ERROR_CONTROL_RESIDUAL = 4;

const double FLOW_TI_ABSOLUTE_TOLERANCE = 1.0; // defaults for time integrations
const double FLOW_TI_RELATIVE_TOLERANCE = 0.0;
const double FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE = 1e-6;
const int FLOW_TI_MAX_ITERATIONS = 400;

const int FLOW_DT_ADAPTIVE = 1;
const double FLOW_DT_ADAPTIVE_INCREASE = 4.0;
const double FLOW_DT_ADAPTIVE_REDUCTION = 0.1;
const double FLOW_DT_ADAPTIVE_SAFETY_FACTOR = 0.9;
const double FLOW_DT_ADAPTIVE_ERROR_TOLERANCE = 1e-10;

// multiscale models
const double FLOW_DPM_NEWTON_TOLERANCE = 1e-8;

const int FLOW_HEX_FACES = 6; // Hexahedron is the common element
const int FLOW_HEX_NODES = 8;
const int FLOW_HEX_EDGES = 12;

const int FLOW_QUAD_FACES = 4; // Quadrilateral is the common element
const int FLOW_QUAD_NODES = 4;
const int FLOW_QUAD_EDGES = 4;

const int FLOW_MAX_FACES = 14; // Kelvin's tetrakaidecahedron
const int FLOW_MAX_NODES = 47; // These polyhedron parameters must
const int FLOW_MAX_EDGES = 60; // be calculated in Init().

const int FLOW_INTERNAL_ERROR = 911; // contact (lipnikov@lanl.gov)

const int FLOW_UPWIND_UPDATE_TIMESTEP = 1;
const int FLOW_UPWIND_UPDATE_ITERATION = 2;

} // namespace Flow
} // namespace Amanzi

#endif
