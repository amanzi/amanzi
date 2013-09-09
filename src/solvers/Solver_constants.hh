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

const int SOLVER_CONVERGENCE_RELATIVE_RHS = 1;  // must be power of 2
const int SOLVER_CONVERGENCE_RELATIVE_RESIDUAL = 2;
const int SOLVER_CONVERGENCE_ABSOLUTE_RESIDUAL = 4;
const int SOLVER_MAKE_ONE_ITERATION = 8;

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif

