/*
  This is the Nonlinear Solver component of the Amanzi code.

  Interface for a nonlinear solver.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#ifndef AMANZI_SOLVER_BASE_
#define AMANZI_SOLVER_BASE_

#include "Teuchos_RCP.hpp"

#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class Solver {
 public:

  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
                    const VectorSpace& map) = 0;

  virtual int Solve(const Teuchos::RCP<Vector>& u) = 0;

  virtual double residual() = 0;
  virtual int num_itrs() = 0;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
