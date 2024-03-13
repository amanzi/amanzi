/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Solvers

*/

#ifndef AMANZI_SOLVER_CONSTANTS_HH_
#define AMANZI_SOLVER_CONSTANTS_HH_

namespace Amanzi {
namespace AmanziSolvers {

// Note: these are powers of 2 to allow multiple convergence criteria, each of
// which can be turned on or off.  Currently convergence is met if any of the
// enabled criteria match.  The exception to this is ONE_ITERATION, which
// requires at least one iteration independent of the initial residual.
const int LIN_SOLVER_RELATIVE_RHS = 1; // must be power of 2
const int LIN_SOLVER_RELATIVE_RESIDUAL = 2;
const int LIN_SOLVER_ABSOLUTE_RESIDUAL = 4;
const int LIN_SOLVER_MAKE_ONE_ITERATION = 8;

const int LIN_SOLVER_NON_SPD_APPLY = -1;
const int LIN_SOLVER_NON_SPD_APPLY_INVERSE = -2;
const int LIN_SOLVER_MAX_ITERATIONS = -3;
const int LIN_SOLVER_RESIDUAL_OVERFLOW = -4;

// Trilinos flags
const int LIN_SOLVER_BELOS_SAYS_SUCCESS = 1;
const int LIN_SOLVER_BELOS_SAYS_FAIL = -1;

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
