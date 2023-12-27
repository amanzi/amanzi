/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Amanzi supports a few nonlinear solvers.
Typically, a process kernel uses a factory to select a nonlinear solver.
This factory uses parameter *solver type* to find parameters for
the selected solver.

*/


#ifndef AMANZI_SOLVER_BASE_
#define AMANZI_SOLVER_BASE_

#include <utility>

#include "Teuchos_RCP.hpp"

#include "ResidualDebugger.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class Solver {
 public:
  virtual ~Solver() = default;

  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn,
                    const Teuchos::RCP<const VectorSpace>& map) = 0;

  virtual int Solve(const Teuchos::RCP<Vector>& u) = 0;

  // mutators
  virtual void set_tolerance(double tol) = 0;
  virtual void set_pc_lag(int pc_lag) = 0;
  virtual void set_db(const Teuchos::RCP<ResidualDebugger>& db) {}

  // access
  virtual double tolerance() = 0;
  virtual double residual() = 0;
  virtual int num_itrs() = 0;
  virtual int returned_code() = 0;
  virtual int pc_calls() = 0;
  virtual int pc_updates() = 0;
};

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
