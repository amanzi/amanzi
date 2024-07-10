/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

*/

#ifndef AMANZI_WHETSTONE_DEFS_HH_
#define AMANZI_WHETSTONE_DEFS_HH_

#include <tuple>
#include <vector>

#include "GeometryDefs.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace WhetStone {

// This definition allows us to use WhetStone as a standalone library.
#define AMANZI_CODE
#define WHETSTONE_ABC 1
#define WHETSTONE_VIEW_SEMANTICS 1

enum class ProjectorType {
  L2,
  H1,
  LS // least square
};

enum class DOF_Type { SCALAR = 1, VECTOR, POINT, NORMAL_COMPONENT, MOMENT };

#ifdef AMANZI_CODE
using Entity_ID = AmanziGeometry::Entity_ID;
using Entity_ID_View = AmanziMesh::Mesh::Entity_ID_View;
using Parallel_kind = AmanziMesh::Parallel_kind;
using Entity_kind = AmanziMesh::Entity_kind;
using SchemaItem = std::tuple<AmanziMesh::Entity_kind, DOF_Type, int>;
using Entity_kind = AmanziMesh::Entity_kind;

#else
typedef long long int Entity_ID;
typedef Entity_ID_List Entity_ID_View;

enum Entity_kind { NODE = 0, EDGE, FACE, CELL, BOUNDARY_FACE };

typedef std::tuple<Entity_kind, DOF_Type, int> SchemaItem;

enum class Parallel_kind {
  OWNED = 1; // Owned by this processor
  GHOST = 2; // Owned by another processor
  ALL = 3;   // OWNED + GHOST
};
#endif

// status of elemental matrices
const int WHETSTONE_ELEMENTAL_MATRIX_OK = 0;
const int WHETSTONE_ELEMENTAL_MATRIX_SIZE = 1;
const int WHETSTONE_ELEMENTAL_MATRIX_FAILED = 2; // only for unexpected situations

// control of the stabilization term in MFD schemes
const int WHETSTONE_STABILITY_GENERIC = 1;
const int WHETSTONE_STABILITY_GENERIC_SCALED = 2;
const int WHETSTONE_STABILITY_OPTIMIZED_DMP = 3;
const int WHETSTONE_STABILITY_OPTIMIZED_GEOMETRY = 4;

const double WHETSTONE_TOLERANCE_DECOMPOSITION = 1e-12;

// control of simplex method
const double WHETSTONE_SIMPLEX_TOLERANCE = 1e-10;
const double WHETSTONE_SIMPLEX_MAX_ITERATIONS = 100; // factor of number of unknowns
const double WHETSTONE_SIMPLEX_NO_FEASIBLE_SET = -1;
const double WHETSTONE_SIMPLEX_NO_CONVERGENCE = -2;
const double WHETSTONE_SIMPLEX_UNBOUNDED_PROBLEM = -3;
const int WHETSTONE_SIMPLEX_FUNCTIONAL_SUMALL = 1;
const int WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE = 2;

#undef WHETSTONE_SIMPLEX_PIVOT_BRANDT // select pivot rule
#define WHETSTONE_SIMPLEX_PIVOT_MFD3D

// various MFD schemes for diffusion
const int DIFFUSION_OPTIMIZED_FOR_SPARSITY = 9; // recommended
const int DIFFUSION_POLYHEDRA_SCALED = 2;
const int DIFFUSION_OPTIMIZED_FOR_MONOTONICITY = 3;
const int DIFFUSION_HEXAHEDRA_MONOTONE = 4;
const int DIFFUSION_SUPPORT_OPERATOR = 7;
const int DIFFUSION_TPFA = 5;

const int ELECTROMAGNETICS_DEFAULT = 1;
const int ELECTROMAGNETICS_GENERALIZED = 2;

// various DG schemes
const int TAYLOR_BASIS_NATURAL = 1;
const int TAYLOR_BASIS_NORMALIZED = 2;
const int TAYLOR_BASIS_NORMALIZED_ORTHO = 3; // recommended
const int TAYLOR_BASIS_REGULARIZED = 4;

} // namespace WhetStone
} // namespace Amanzi

#endif
