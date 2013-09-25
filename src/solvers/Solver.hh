/*
  This is the Nonlinear Solver component of the Amanzi code.

  Interface for a nonlinear solver.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#ifndef AMANZI_SOLVER_BASE_
#define AMANZI_SOLVER_BASE_

#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector>
class Solver {
 public:
  virtual int Solve(const Teuchos::RCP<Vector>& u) = 0;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
