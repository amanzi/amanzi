/*
This is the Linear Solver component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_SOLVER_CONSTANTS_HH_
#define AMANZI_SOLVER_CONSTANTS_HH_

namespace Amanzi {
namespace AmanziSolvers {

const int LIN_SOLVER_RELATIVE_RHS = 1;  // must be power of 2
const int LIN_SOLVER_RELATIVE_RESIDUAL = 2;
const int LIN_SOLVER_ABSOLUTE_RESIDUAL = 4;
const int LIN_SOLVER_MAKE_ONE_ITERATION = 8;

const int LIN_SOLVER_NON_SPD_APPLY = -1;
const int LIN_SOLVER_NON_SPD_APPLY_INVERSE = -1;
const int LIN_SOLVER_MAX_ITERATIONS = -2;
const int LIN_SOLVER_RESIDUAL_OVERFLOW = -3;

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif

