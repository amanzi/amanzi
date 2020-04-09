/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Interface for a Nonlinear Solver

#ifndef AMANZI_SOLVER_
#define AMANZI_SOLVER_

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

  // Returns 0 if success, 1 if failure.
  virtual int Solve(const Teuchos::RCP<Vector>& u) = 0;

  // mutators
  virtual void set_tolerance(double tol) = 0;
  virtual void set_pc_lag(int pc_lag) = 0;
  virtual void set_db(const Teuchos::RCP<ResidualDebugger>& db) = 0;

  // accessors

  // Note, what the error is is dependent upon the MonitorType and the
  // MonitorNorm
  virtual double error() const = 0;

  // Tolerance to compare to error.
  virtual double tolerance() const = 0;

  // The L2 norm of the residual
  virtual double residual() const = 0;

  // Number of nonlinear iterations
  virtual int num_iterations() const = 0;

  // See SolverDefs.hh MonitorStatus definition, but positive values indicate
  // number of iterations convergence was achieved in while negative numbers
  // indicate an error.
  virtual int returned_code() const = 0;

  // Number of preconditioner ApplyInverse calls (this solve)
  virtual int pc_calls() const = 0;

  // Number of times predonditioner was updated (this solve)
  virtual int pc_updates() const = 0;

  // Number of times nonlinear residual function was called (this solve)
  virtual int function_calls() const = 0;

  // name of the method (for logging)
  virtual std::string name() const = 0;
};

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
