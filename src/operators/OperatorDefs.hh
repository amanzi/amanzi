/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATORS_DEFS_HH_
#define AMANZI_OPERATORS_DEFS_HH_

namespace Amanzi {
namespace Operators {

typedef enum {
  OPERATOR_DIFFUSION_MFD,
  OPERATOR_DIFFUSION_FV,
  OPERATOR_DIFFUSION_NLFV,
  OPERATOR_DIFFUSION_NLFVFACES,
  OPERATOR_DIFFUSION_MFD_GRAVITY,
  OPERATOR_DIFFUSION_FV_GRAVITY,
  OPERATOR_DIFFUSION_NLFV_GRAVITY,
  OPERATOR_DIFFUSION_NLFVFACES_GRAVITY,
  OPERATOR_DIFFUSION_MFD_XMOF,
  OPERATOR_DIFFUSION_MFD_TRACER,
  OPERATOR_DIFFUSION_DG,
  OPERATOR_DIFFUSION_FRACTURED_MATRIX,
  OPERATOR_ADVECTION,
  OPERATOR_ACCUMULATION,
  OPERATOR_ELASTICITY,
  OPERATOR_ELECTROMAGNETICS,
  OPERATOR_MAGNETIC_DIFFUSION
} OperatorType;

enum class DOF_Type { SCALAR = 1, VECTOR, POINT, NORMAL_COMPONENT, MOMENT };

// Constants in the next block must powers of 2.
const int OPERATOR_SCHEMA_DOFS_FACE = 1;
const int OPERATOR_SCHEMA_DOFS_CELL = 2;
const int OPERATOR_SCHEMA_DOFS_NODE = 4;
const int OPERATOR_SCHEMA_DOFS_EDGE = 8;
const int OPERATOR_SCHEMA_DOFS_BNDFACE = 16;

const int OPERATOR_SCHEMA_BASE_FACE = 32;
const int OPERATOR_SCHEMA_BASE_CELL = 64;
const int OPERATOR_SCHEMA_BASE_NODE = 128;
const int OPERATOR_SCHEMA_BASE_EDGE = 256;

const int OPERATOR_SCHEMA_INDICES = 512;

// schemas
const int OPERATOR_SCHEMA_RULE_EXACT = 1;
const int OPERATOR_SCHEMA_RULE_SUBSET = 2;

// Boundary Conditions:
//   Dirichlet, Neumann and Mixed are conventional boundary conditions
//   for 2nd-order operators. Composite (additive) operators may require
//   special treatment of total flux conditions. Finally some essential
//   boundary conditions may be imposed in a weak form which leads to
//   type2 boundary conditions. See BCs.hh for more detail.
const int OPERATOR_BC_NONE = 0;
const int OPERATOR_BC_DIRICHLET = 1;
const int OPERATOR_BC_NEUMANN = 2;
const int OPERATOR_BC_TOTAL_FLUX = 3;
const int OPERATOR_BC_MIXED = 4;
const int OPERATOR_BC_DIRICHLET_TYPE2 = 5;
const int OPERATOR_BC_REMOVE = 6;

// memory allocation
const int OPERATOR_HEX_FACES = 6; // Hexahedron is the common element
const int OPERATOR_HEX_NODES = 8;
const int OPERATOR_HEX_EDGES = 12;

const int OPERATOR_QUAD_FACES = 4; // Quadrilateral is the common element
const int OPERATOR_QUAD_NODES = 4;
const int OPERATOR_QUAD_EDGES = 4;

// Newton-correction options
const int OPERATOR_DIFFUSION_JACOBIAN_NONE = 0;
const int OPERATOR_DIFFUSION_JACOBIAN_TRUE = 1;
const int OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE = 2;

// upwind options
const int OPERATOR_UPWIND_NONE = 0;
const int OPERATOR_UPWIND_CONSTANT_VECTOR = 1;
const int OPERATOR_UPWIND_FLUX = 2;
const int OPERATOR_UPWIND_GRAVITY = 4;
const int OPERATOR_UPWIND_DIVK = 8;
const int OPERATOR_UPWIND_ARITHMETIC_AVERAGE = 16;
const int OPERATOR_UPWIND_FLUX_SECOND_ORDER = 32;
const double OPERATOR_UPWIND_RELATIVE_TOLERANCE = 1e-12;

// method for nonlinear coefficient (use power of 2)
const int OPERATOR_LITTLE_K_NONE = 0;
const int OPERATOR_LITTLE_K_UPWIND = 1;
const int OPERATOR_LITTLE_K_DIVK_BASE = 2;       // base (only face component)
const int OPERATOR_LITTLE_K_DIVK = 6;            // add cell component
const int OPERATOR_LITTLE_K_DIVK_TWIN = 10;      // add twin component
const int OPERATOR_LITTLE_K_DIVK_TWIN_GRAD = 18; // add grad component
const int OPERATOR_LITTLE_K_STANDARD = 32;

// method for gravity
const int OPERATOR_GRAVITY_HH = 1;
const int OPERATOR_GRAVITY_FV = 2;

// special properties of operators
const int OPERATOR_PROPERTY_DATA_READ_ONLY = 1; // must be power of 2
const int OPERATOR_PROPERTY_MAP_CONSTANT = 2;

// reconstruction options
const double OPERATOR_RECONSTRUCTION_MATRIX_CORRECTION = 1e-15;

// limiting options
const int OPERATOR_LIMITER_BARTH_JESPERSEN = 1;
const int OPERATOR_LIMITER_BARTH_JESPERSEN_DG = 2;
const int OPERATOR_LIMITER_MICHALAK_GOOCH = 3;
const int OPERATOR_LIMITER_MICHALAK_GOOCH_DG = 4;
const int OPERATOR_LIMITER_TENSORIAL = 5;
const int OPERATOR_LIMITER_KUZMIN = 6;
const double OPERATOR_LIMITER_TOLERANCE = 1e-14;
const double OPERATOR_LIMITER_INFINITY = 1e+99;

// stencil for calculating limiting bounds
const int OPERATOR_LIMITER_STENCIL_N2C = 10;
const int OPERATOR_LIMITER_STENCIL_F2C = 20;
const int OPERATOR_LIMITER_STENCIL_C2C_CLOSEST = 30;
const int OPERATOR_LIMITER_STENCIL_C2C_ALL = 31;

// detrminant of Jacobian of a map
const int OPERATOR_DETERMINANT_VEM = 1;
const int OPERATOR_DETERMINANT_EXACT_TI = 2; // experimental
const int OPERATOR_DETERMINANT_MONOTONE = 3;

} // namespace Operators
} // namespace Amanzi

#endif
