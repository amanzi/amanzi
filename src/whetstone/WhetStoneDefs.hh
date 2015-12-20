/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_WHETSTONE_DEFS_HH_
#define AMANZI_WHETSTONE_DEFS_HH_

namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_ELEMENTAL_MATRIX_OK = 0;
const int WHETSTONE_ELEMENTAL_MATRIX_WRONG = 1;
const int WHETSTONE_ELEMENTAL_MATRIX_PASSED = 2;
const int WHETSTONE_ELEMENTAL_MATRIX_FAILED = 4;  // only for unexpected situations

const int WHETSTONE_STABILITY_GENERIC = 1;
const int WHETSTONE_STABILITY_GENERIC_SCALED = 2;
const int WHETSTONE_STABILITY_OPTIMIZED_DMP = 3;
const int WHETSTONE_STABILITY_OPTIMIZED_GEOMETRY = 4;

const int WHETSTONE_MAX_SPATIAL_DIMENSION = 3;

const double WHETSTONE_TOLERANCE_DECOMPOSITION = 1e-12;

const double WHETSTONE_SIMPLEX_TOLERANCE = 1e-10;
const double WHETSTONE_SIMPLEX_MAX_ITERATIONS = 100;  // factor of number of unknowns
const double WHETSTONE_SIMPLEX_NO_FEASIBLE_SET = -1;
const double WHETSTONE_SIMPLEX_NO_CONVERGENCE = -2;
const double WHETSTONE_SIMPLEX_UNBOUNDED_PROBLEM = -3;
const int WHETSTONE_SIMPLEX_FUNCTIONAL_SUMALL = 1;
const int WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE = 2;

#undef WHETSTONE_SIMPLEX_PIVOT_BRANDT  // select pivot rule
#define WHETSTONE_SIMPLEX_PIVOT_MFD3D

const int DIFFUSION_OPTIMIZED_FOR_SPARSITY = 9;  // recommended MFD methods
const int DIFFUSION_POLYHEDRA_SCALED = 2; 
const int DIFFUSION_OPTIMIZED_FOR_MONOTONICITY = 3;  
const int DIFFUSION_HEXAHEDRA_MONOTONE = 4;
const int DIFFUSION_SUPPORT_OPERATOR = 7;
const int DIFFUSION_TPFA = 5; 

}  // namespace WhetStone
}  // namespace Amanzi

#endif

