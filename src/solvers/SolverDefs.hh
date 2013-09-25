/*
This is the Nonlinear Solver component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_SOLVER_DEFS_HH_
#define AMANZI_SOLVER_DEFS_HH_

namespace Amanzi {
namespace AmanziSolvers {

enum ConvergenceMonitor {
     SOLVER_MONITOR_UPDATE = 0,
     SOLVER_MONITOR_PCED_RESIDUAL = 1,
     SOLVER_MONITOR_RESIDUAL = 2
};

const int SOLVER_CONTINUE = 1;
const int SOLVER_CONVERGED = 0;

const int SOLVER_MAX_ITERATIONS = -1;
const int SOLVER_OVERFLOW = -2;
const int SOLVER_STAGNATING = -3;
const int SOLVER_DIVERGING = -4;

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif

