/*
This is the Nonlinear Solver component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_NONLINEAR_SOLVER_CONSTANTS_HH_
#define AMANZI_NONLINEAR_SOLVER_CONSTANTS_HH_

namespace Amanzi {
namespace AmanziSolvers {

const int SOLVER_CONVERGED = 0;

const int SOLVER_MAX_ITERATIONS = -1;
const int SOLVER_OVERFLOW = -2;
const int SOLVER_STAGNATING = -3;
const int SOLVER_DIVERGING = -4;

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif

