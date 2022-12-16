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

  Interface for a nonlinear solver.
*/


#ifndef AMANZI_SOLVER_BASE_
#define AMANZI_SOLVER_BASE_

#include "Teuchos_RCP.hpp"

#include "ResidualDebugger.hh"
#include "SolverDefs.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class Solver {
 public:
  virtual ~Solver() = default;

  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn, const VectorSpace& map) = 0;

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


// non-member functions for parsing input plist
inline void
ParseConvergenceCriteria(const std::string& monitor_name,
                         ConvergenceMonitor* monitor,
                         int* norm_type)
{
  *norm_type = SOLVER_NORM_LINF;
  if (monitor_name == "monitor residual") {
    *monitor = SOLVER_MONITOR_RESIDUAL;
  } else if (monitor_name == "monitor l2 residual") {
    *monitor = SOLVER_MONITOR_RESIDUAL;
    *norm_type = SOLVER_NORM_L2;
  } else if (monitor_name == "monitor preconditioned residual") {
    *monitor = SOLVER_MONITOR_PCED_RESIDUAL;
  } else if (monitor_name == "monitor preconditioned l2 residual") {
    *monitor = SOLVER_MONITOR_PCED_RESIDUAL;
    *norm_type = SOLVER_NORM_L2;
  } else if (monitor_name == "monitor update") {
    *monitor = SOLVER_MONITOR_UPDATE; // default value
  } else {
    Errors::Message m;
    m << "Invalid monitor name for nonlinear solver: \"" << monitor_name << "\"";
    Exceptions::amanzi_throw(m);
  }
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
