/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __FLOW_CONSTANTS_HH__
#define __FLOW_CONSTANTS_HH__


namespace Amanzi {
namespace AmanziFlow {

const int FLOW_STATUS_NULL = 0;  // used for internal debuging
const int FLOW_STATUS_INIT = 2;
const int FLOW_STATUS_INITIAL_GUESS = 4;
const int FLOW_STATUS_STEADY_STATE = 4;
const int FLOW_STATUS_TRANSIENT_STATE = 8;

const int FLOW_STATE_VIEW = 1;  // copy Teuchos::RCP pointers
const int FLOW_STATE_COPY = 2;  // add ghost data to some arrays

const int FLOW_BC_FACE_NULL = 0; 
const int FLOW_BC_FACE_PRESSURE = 1; 
const int FLOW_BC_FACE_PRESSURE_SEEPAGE = 2; 
const int FLOW_BC_FACE_FLUX = 3;
const int FLOW_BC_FACE_MIXED = 4;

const int FLOW_BC_SUBMODEL_RAINFALL = 1;
const int FLOW_BC_SUBMODEL_SEEPAGE_PFLOTRAN = 2;
const int FLOW_BC_SUBMODEL_SEEPAGE_FACT = 4;
const int FLOW_BC_SUBMODEL_HEAD_RELATIVE = 8;
const double FLOW_BC_SEEPAGE_FACE_REGULARIZATION = 100;  // [Pa]

const int FLOW_SOLVER_NKA = 1;
const int FLOW_SOLVER_NEWTON = 2;
const int FLOW_SOLVER_PICARD_NEWTON = 3;

const int FLOW_TIME_INTEGRATION_PICARD = 1;
const int FLOW_TIME_INTEGRATION_BACKWARD_EULER = 2;  // Only for testing.
const int FLOW_TIME_INTEGRATION_BDF1 = 3;
const int FLOW_TIME_INTEGRATION_BDF2 = 4;
const double FLOW_INITIAL_DT = 1e-8;  // [sec]
const double FLOW_MAXIMUM_DT = 3.15e+10;  // [sec] 1000 years
const double FLOW_YEAR = 3.15576e+7;

const int FLOW_RELATIVE_PERM_NONE = 1; 
const int FLOW_RELATIVE_PERM_CENTERED = 2; 
const int FLOW_RELATIVE_PERM_UPWIND_GRAVITY = 3; 
const int FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX = 4;
const int FLOW_RELATIVE_PERM_ARITHMETIC_MEAN = 5;
const int FLOW_RELATIVE_PERM_EXPERIMENTAL = 6;
const double FLOW_RELATIVE_PERM_TOLERANCE = 1e-10;  // [-]

const int FLOW_MFD3D_POLYHEDRA = 1;  // default
const int FLOW_MFD3D_HEXAHEDRA_MONOTONE = 3;  // highly experimental
const int FLOW_MFD3D_TWO_POINT_FLUX = 4;  // without consistency
const int FLOW_MFD3D_SUPPORT_OPERATOR = 5;
const int FLOW_MFD3D_OPTIMIZED = 6;
const int FLOW_MFD3D_OPTIMIZED_EXPERIMENTAL = 7;  // experimental

const int FLOW_PRECONDITIONER_TRILINOS_ML = 1;  // preconditioners
const int FLOW_PRECONDITIONER_HYPRE_AMG = 2;
const int FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU = 3;

const int FLOW_TI_ERROR_CONTROL_PRESSURE = 1;  // binary mask for error control
const int FLOW_TI_ERROR_CONTROL_SATURATION = 2;
const int FLOW_TI_ERROR_CONTROL_RESIDUAL = 4;

const double FLOW_TI_ABSOLUTE_TOLERANCE = 1.0;  // defaults for time integrations
const double FLOW_TI_RELATIVE_TOLERANCE = 0.0;
const double FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE = 1e-6;
const int FLOW_TI_MAX_ITERATIONS = 400;

const int FLOW_MATRIX_MFD = 1;  // matrix to use in linear and nonlinear solvers 
const int FLOW_MATRIX_MFD_TPFA = 2;

const int FLOW_MATRIX_ACTION_MATRIX = 1;
const int FLOW_MATRIX_ACTION_PRECONDITIONER = 2;

const int FLOW_DT_ADAPTIVE = 1;
const double FLOW_DT_ADAPTIVE_INCREASE = 4.0;
const double FLOW_DT_ADAPTIVE_REDUCTION = 0.1;
const double FLOW_DT_ADAPTIVE_SAFETY_FACTOR = 0.9;
const double FLOW_DT_ADAPTIVE_ERROR_TOLERANCE = 1e-10;

const int FLOW_HEX_FACES = 6;  // Hexahedron is the common element
const int FLOW_HEX_NODES = 8;
const int FLOW_HEX_EDGES = 12;

const int FLOW_QUAD_FACES = 4;  // Quadrilateral is the common element
const int FLOW_QUAD_NODES = 4;
const int FLOW_QUAD_EDGES = 4;

const int FLOW_MAX_FACES = 14;  // Kelvin's tetrakaidecahedron
const int FLOW_MAX_NODES = 47;  // These polyhedron parameters must
const int FLOW_MAX_EDGES = 60;  // be calculated in Init().

const int FLOW_INTERNAL_ERROR = 911;  // contact (lipnikov@lanl.gov)

const int FLOW_VERBOSITY_NONE = 0;
const int FLOW_VERBOSITY_LOW = 1;
const int FLOW_VERBOSITY_MEDIUM = 2;
const int FLOW_VERBOSITY_HIGH = 3;
const int FLOW_VERBOSITY_EXTREME = 4;

const int FLOW_AMANZI_VERSION = 2;  

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

