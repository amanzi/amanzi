/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATORS_DEFS_HH_
#define AMANZI_OPERATORS_DEFS_HH_

namespace Amanzi {
namespace Operators {

// general information about an operator, e.g. a preconditioner may
// be wrapped up in an iterative solver or be of a "raw" matrix type
typedef enum { OPERATOR_MATRIX,
               OPERATOR_PRECONDITIONER,
               OPERATOR_PRECONDITIONER_RAW,
               OPERATOR_TERM_DIFFUSION } OperatorType;

// this is not used currently and my go away
typedef enum { PDE_DIFFUSION,
               PDE_DIFFUSION_MFD,
               PDE_DIFFUSION_FV,
               PDE_DIFFUSION_NLFV,
               PDE_DIFFUSION_NLFVFACES,
               PDE_DIFFUSION_MFD_GRAVITY,
               PDE_DIFFUSION_FV_GRAVITY,
               PDE_DIFFUSION_NLFV_GRAVITY,
               PDE_DIFFUSION_NLFVFACES_GRAVITY,
               PDE_DIFFUSION_MFD_XMOF,
               PDE_DIFFUSION_MFD_TRACER,
               PDE_DIFFUSION_DG,
               PDE_DIFFUSION_FRACTURED_MATRIX,
               PDE_ADVECTION,
               PDE_ACCUMULATION,
               PDE_ELASTICITY,
               PDE_ELECTROMAGNETICS,
               PDE_MAGNETIC_DIFFUSION} PDEType;

// coefficient type
typedef enum { CONSTANT = 0,  // includes tensorial coefficients
               POLYNOMIAL,
               VECTOR_POLYNOMIAL,
               VECTOR_SPACETIME_POLYNOMIAL,
               MATRIX_POLYNOMIAL,
               FUNCTION } CoefType;

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
const int OPERATOR_SCHEMA_RULE_SUPERSET = 3;

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
const int OPERATOR_HEX_FACES = 6;  // Hexahedron is the common element
const int OPERATOR_HEX_NODES = 8;
const int OPERATOR_HEX_EDGES = 12;

const int OPERATOR_QUAD_FACES = 4;  // Quadrilateral is the common element
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
const int OPERATOR_LITTLE_K_DIVK_BASE = 2;  // base (only face component)
const int OPERATOR_LITTLE_K_DIVK = 6;  // add cell component
const int OPERATOR_LITTLE_K_DIVK_TWIN = 10;  // add twin component
const int OPERATOR_LITTLE_K_STANDARD = 32;

// method for gravity
const int OPERATOR_GRAVITY_HH = 1;
const int OPERATOR_GRAVITY_FV = 2;

// reconstruction options
const double OPERATOR_RECONSTRUCTION_MATRIX_CORRECTION = 1e-15;

// limiting options
const int OPERATOR_LIMITER_BARTH_JESPERSEN = 1;
const int OPERATOR_LIMITER_BARTH_JESPERSEN_DG = 2;
const int OPERATOR_LIMITER_MICHALAK_GOOCH = 3;
const int OPERATOR_LIMITER_MICHALAK_GOOCH_DG = 4;
const int OPERATOR_LIMITER_TENSORIAL = 5;
const int OPERATOR_LIMITER_KUZMIN = 6;
const int OPERATOR_LIMITER_BARTH_JESPERSEN_DG_HIERARCHICAL = 7;

const double OPERATOR_LIMITER_TOLERANCE = 1e-14;
const double OPERATOR_LIMITER_INFINITY = 1e+99;

// stencil for calculating limiting bounds
const int OPERATOR_LIMITER_STENCIL_N2C = 10;
const int OPERATOR_LIMITER_STENCIL_E2C = 20;
const int OPERATOR_LIMITER_STENCIL_F2C = 30;
const int OPERATOR_LIMITER_STENCIL_C2C_CLOSEST = 40;
const int OPERATOR_LIMITER_STENCIL_C2C_ALL = 41;

const int OPERATOR_MAX_NUM_FACES = 10; // wild guess at this point...

}  // namespace Operators
}  // namespace Amanzi

#endif

