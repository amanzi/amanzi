/*
  This is the Nonlinear Solver component of the Amanzi code.

  Interface for a nonlinear solver.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#ifndef AMANZI_SOLVER_BASE_
#define AMANZI_SOLVER_BASE_

#include "Teuchos_RCP.hpp"

#include "ResidualDebugger.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class Solver {
 public:

  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
                    const VectorSpace& map) = 0;

  virtual int Solve(const Teuchos::RCP<Vector>& u) = 0;

  // mutators
  virtual void set_tolerance(double tol) = 0;
  virtual void set_pc_lag(double pc_lag) = 0;
  virtual void set_db(const Teuchos::RCP<
		      ResidualDebugger<Vector,VectorSpace> >& db) {}
  
  // access 
  virtual double tolerance() = 0;
  virtual double residual() = 0;
  virtual int num_itrs() = 0;
  virtual int returned_code() = 0;
  virtual int pc_calls() = 0;
  virtual int pc_updates() = 0;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
