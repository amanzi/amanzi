/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Solves globally implicit systems using backward Euler
/*!

Backward Euler is the simplest of the implicit methods.  It solves time
integration schemes by evaluating all time derivatives at the new time.  This
makes it unconditionally stable, though potentially not very accurate.  This
unconditional stability tends to make it the workhorse of the types of stiff,
nonlinear parabolic equations such as Richards equation and the diffusion wave
approximation.

In this method, we look to solve:

.. math::
    \frac{\partial \mathbf{u}}{\partial t} = f(\mathbf{u},\mathbf{x},t)

via the time discretization scheme:

.. math::
    \frac{\mathbf{u}^{t + \Delta t} - \mathbf{u}^{t}}{\Delta t} = f(\mathbf{u}^{t + \Delta t}, \mathbf{x}, t + \Delta t)

.. _bdf1-ti-spec:
.. admonition:: bdf1-ti-spec

   * `"verbose object`" ``[verbose-object-spec]`` A `Verbose Object`_

   * `"residual debugger`" ``[residual-debugger-spec]`` A `Residual Debugger`_ object.

   * `"max preconditioner lag iterations`" ``[int]`` **0** specifies frequency
     of preconditioner recalculation.

   * `"freeze preconditioner`" ``[bool]`` **false** enforces preconditioner to
     be updated only once per non-linear solver. When set to true, the above
     parameter is ignored.

   * `"extrapolate initial guess`" ``[bool]`` **true** identifies forward time
     extrapolation of the initial guess.

   * `"nonlinear iteration initial guess extrapolation order`" ``[int]`` **1**
     defines extrapolation algorithm. Zero value implies no extrapolation.

   * `"restart tolerance relaxation factor`" ``[double]`` **1** Changes the
     nonlinear tolerance on restart. The time integrator is usually restarted
     when a boundary condition changes drastically. It may be beneficial to
     loosen the nonlinear tolerance on the first several timesteps after the
     time integrator restart. The default value is 1, while a reasonable value
     may be as large as 1000.

   * `"restart tolerance relaxation factor damping`" ``[double]`` **1**
     Controls how fast the loosened nonlinear tolerance will revert back to the
     one specified in `"nonlinear tolerance`". If the nonlinear tolerance is
     "tol", the relaxation factor is "factor", and the damping is "d", and the
     timestep count is "n" then the actual nonlinear tolerance is "tol *
     max(1.0, factor * d ** n)". Reasonable values are between 0 and 1.

   INCLUDES

   - ``[solver-typed-spec]`` *Uses a* Solver_.
   - ``[timestep-controller-typed-spec]`` *Uses a* `Timestep Controller`_


Note this also accepts an object that provides the `BDF1 Solver Interface`_.

.. code-block:: xml

  <ParameterList name="time integrator">
    <Parameter name="time integration method" type="string" value="BDF1"/>
    <ParameterList name="BDF1">
      <Parameter name="max preconditioner lag iterations" type="int" value="5"/>
      <Parameter name="freeze preconditioner" type="bool" value="false"/>
      <Parameter name="extrapolate initial guess" type="bool" value="true"/>
      <Parameter name="nonlinear iteration initial guess extrapolation order" type="int" value="1"/>
      <Parameter name="restart tolerance relaxation factor" type="double" value="1.0"/>
      <Parameter name="restart tolerance relaxation factor damping" type="double" value="1.0"/>

      <Parameter name="timestep controller type" type="string" value="standard"/>
      <ParameterList name="timestep controller standard parameters">
        ...
      </ParameterList>

      <Parameter name="solver type" type="string" value="nka"/>
      <ParameterList name="nka parameters">
        ...
      </ParameterList>
    </ParameterList>
  </ParameterList>

*/


#ifndef AMANZI_BDF1_TIME_INTEGRATOR_HH_
#define AMANZI_BDF1_TIME_INTEGRATOR_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "errors.hh"
#include "VerboseObject.hh"
#include "Solver.hh"
#include "SolverFactory.hh"
#include "SolverDefs.hh"

#include "BDF1_State.hh"
#include "BDF1_SolverFnBase.hh"

namespace Amanzi {

template<class Vector, class VectorSpace>
class BDF1_TI {
 public:
  // Create the BDF Dae solver object, the nonlinear problem must be
  // defined in a class that derives from the virtual base class FnBase.
  BDF1_TI(const std::string& name,
          Teuchos::ParameterList& plist,
          BDFFnBase<Vector>& fn,
          const Teuchos::RCP<const VectorSpace>& space,
          const Teuchos::RCP<State>& S = Teuchos::null);

  // initializes the state
  void SetInitialState(const double h,
                       const Teuchos::RCP<const Vector>& u,
                       const Teuchos::RCP<const Vector>& udot = Teuchos::null);

  // After a successful step, this method commits the new
  // solution to the solution history
  void CommitSolution(const double h, const Teuchos::RCP<const Vector>& u);

  // Computes a step and returns true whan it fails.
  bool AdvanceStep(double dt,
                   const Teuchos::RCP<const Vector>& u_prev,
                   const Teuchos::RCP<Vector>& u,
                   double& dt_next);

  bool AdvanceStep(double dt, double& dt_next, const Teuchos::RCP<Vector>& x)
  {
    return AdvanceStep(dt, state_->uhist->MostRecentSolution(), x, dt_next);
  }

  // Reset the memory of the time integrator
  void Reset();

  // returns the most recent time
  double time();

  // returns the initial step size
  double initial_timestep() { return ts_control_->getInitialTimestep(); }

  // returns current nonlinear tolerance
  double tol_solver() { return tol_solver_; }

  // Report statistics
  int number_nonlinear_steps() { return state_->solve_itrs; }
  void ReportStatistics_();
  int number_solver_iterations() { return solver_->num_itrs(); }

 protected:
  Teuchos::RCP<TimestepController> ts_control_; // timestep controller
  Teuchos::RCP<BDF1_State<Vector, VectorSpace>> state_;

  Teuchos::RCP<AmanziSolvers::Solver<Vector, VectorSpace>> solver_;
  Teuchos::RCP<BDF1_SolverFnBase<Vector>> solver_fn_;
  Teuchos::RCP<BDFFnBase<Vector>> fn_;

  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<AmanziSolvers::ResidualDebugger> db_;

  Teuchos::RCP<Vector> udot_prev_, udot_; // for error estimate

 private:
  double tol_solver_; // reference solver's tolerance
};


/* ******************************************************************
* Constructor
****************************************************************** */
template<class Vector, class VectorSpace>
BDF1_TI<Vector, VectorSpace>::BDF1_TI(const std::string& name,
                                      Teuchos::ParameterList& plist,
                                      BDFFnBase<Vector>& fn,
                                      const Teuchos::RCP<const VectorSpace>& space,
                                      const Teuchos::RCP<State>& S)
{
  fn_ = Teuchos::rcpFromRef(fn);

  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject(space->Comm(), name, plist));
  db_ = Teuchos::rcp(new AmanziSolvers::ResidualDebugger(plist.sublist("residual debugger"), S));

  // Create the state.
  state_ = Teuchos::rcp(new BDF1_State<Vector, VectorSpace>(name, plist, space, S));

  // Set up the nonlinear solver
  // -- initialized the SolverFnBase interface
  solver_fn_ = Teuchos::rcp(new BDF1_SolverFnBase<Vector>(plist, fn_));

  AmanziSolvers::SolverFactory<Vector, VectorSpace> factory;
  solver_ = factory.Create(plist);
  solver_->set_db(db_);
  solver_->Init(solver_fn_, *space);

  // Allocate memory for adaptive timestep controll
  udot_ = Teuchos::rcp(new Vector(*space));
  udot_prev_ = Teuchos::rcp(new Vector(*space));

  // timestep controller (note, pointer copy)
  Teuchos::RCP<const Vector> udot_c(udot_);
  Teuchos::RCP<const Vector> udot_prev_c(udot_prev_);
  ts_control_ = createTimestepController(name, plist, S, udot_c, udot_prev_c);

  // misc internal parameters
  tol_solver_ = solver_->tolerance();
}


/* ******************************************************************
* Initialize miscaleneous parameters.
****************************************************************** */
template<class Vector, class VectorSpace>
void
BDF1_TI<Vector, VectorSpace>::SetInitialState(const double t,
                                              const Teuchos::RCP<const Vector>& x,
                                              const Teuchos::RCP<const Vector>& xdot)
{
  // set a clean initial state for when the time integrator is reinitialized
  state_->uhist->FlushHistory(t, *x, xdot.get());
  state_->seq = 0;
  state_->pc_lag = 0;
}


/* ******************************************************************
* Record a successful solution to the history.
****************************************************************** */
template<class Vector, class VectorSpace>
void
BDF1_TI<Vector, VectorSpace>::CommitSolution(const double h, const Teuchos::RCP<const Vector>& u)
{
  double t = h + state_->uhist->MostRecentTime();

  // record the solution for later use when computing an initial guess
  // for the nonlinear solver
  state_->uhist->RecordSolution(t, *u);

  // record some information about this timestep
  state_->hlast = h;
  state_->seq++;
  state_->failed_current = 0;
  state_->hmin = std::min<double>(h, state_->hmin);
  state_->hmax = std::max<double>(h, state_->hmax);
}


/* ******************************************************************
* Returns most recent time.
****************************************************************** */
template<class Vector, class VectorSpace>
double
BDF1_TI<Vector, VectorSpace>::time()
{
  return state_->uhist->MostRecentTime();
}


/* ******************************************************************
* Implementation of implicit Euler timestep.
****************************************************************** */
template<class Vector, class VectorSpace>
bool
BDF1_TI<Vector, VectorSpace>::AdvanceStep(double dt,
                                          const Teuchos::RCP<const Vector>& u_prev,
                                          const Teuchos::RCP<Vector>& u,
                                          double& dt_next)
{
  // initialize the output stream
  Teuchos::OSTab tab = vo_->getOSTab();

  // print info about the timestep
  double tlast = state_->uhist->MostRecentTime();
  double tnew = tlast + dt;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "step " << state_->seq << " T = " << tlast << " [sec]  dT = " << dt << std::endl;
    *vo_->os() << "preconditioner lag is " << state_->pc_lag << " out of " << state_->maxpclag
               << std::endl;
  }

  // Predicted solution (initial value for the nonlinear solver)
  if (state_->extrapolate_guess) {
    if (state_->uhist->history_size() > 1) {
      state_->uhist->InterpolateSolution(tnew, *u);
      fn_->ChangedSolution();

      if (fn_->IsAdmissible(u)) {
        bool changed = fn_->ModifyPredictor(dt, u_prev, u);
        if (changed) fn_->ChangedSolution();
      } else {
        *u = *u_prev;
        fn_->ChangedSolution();
      }
    }
  }

  // Set up the solver fn.
  solver_fn_->SetTimes(tlast, tnew);
  solver_fn_->SetPreviousTimeSolution(u_prev);

  // Set up tolerance due to damping.
  double factor = state_->tol_multiplier;
  double tol = tol_solver_ * factor;
  solver_->set_tolerance(tol);

  if (factor > 1.0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "modified tolerance=" << tol << ", setting PC lag to 0" << std::endl;
    }
  }

  // Update the debugger
  db_->StartIteration<VectorSpace>(state_->failed_current, u->Map());

  // Overwrite preconditioner control
  if (state_->freeze_pc) solver_->set_pc_lag(1000000000);

  // Solve the nonlinear BCE system.
  int ierr, code, itr;
  try {
    ierr = solver_->Solve(u);
    itr = solver_->num_itrs();
    code = solver_->returned_code();
  } catch (const Errors::CutTimestep& e) {
    ierr = 1;
    itr = -1; // This should not be summed up into the global counter.
    code = AmanziSolvers::SOLVER_INTERNAL_EXCEPTION;
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << e.what() << std::endl;
    }
  }

  if (ierr == 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "success: " << itr << " nonlinear itrs"
                 << " error=" << solver_->residual() << std::endl;
    }
  } else {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << vo_->color("red") << "step failed with error code " << code << vo_->reset()
                 << std::endl;
    }
    *u = *u_prev;
  }

  // update the next timestep size
  if (ierr != 0) itr = -1;
  bool is_valid = fn_->IsValid(u);
  dt_next = ts_control_->getTimestep(dt, itr, is_valid);

  // update the preconditioner lag and tolerance multiplier
  if (ierr != 0) {
    state_->pc_lag = 0;
  } else {
    if (dt_next > dt) {
      state_->pc_lag = std::min(state_->pc_lag + 1, state_->maxpclag);
      state_->tol_multiplier = std::max(1.0, factor * state_->tol_multiplier_damp);
    } else if (dt_next < dt) {
      state_->pc_lag = std::max(state_->pc_lag - 1, 0);
    }
  }
  if (factor > 1.0) state_->pc_lag = 0;
  solver_->set_pc_lag(state_->pc_lag);

  // update performance statistics
  state_->solve_itrs += solver_->num_itrs();
  state_->pc_updates += solver_->pc_updates();
  state_->pc_calls += solver_->pc_calls();

  if (ierr != 0) {
    state_->failed_solve++;
    state_->failed_current++;
  } else {
    state_->hmax = std::max(state_->hmax, dt);
    state_->hmin = std::min(state_->hmin, dt);

    if (state_->uhist->history_size() > 1) {
      *udot_prev_ = *udot_;
      double tmp = 1.0 / dt;
      *udot_ = *u;
      udot_->Update(-tmp, *u_prev, tmp);
    }

    ReportStatistics_();
  }

  // debug tool: forcing step to report failure
  if (state_->report_failure == state_->seq) {
    ierr = 1;
    state_->report_failure = -1;
  }

  bool failed = (ierr != 0) || !is_valid;
  return failed;
}


/* ******************************************************************
* Report statistics.
****************************************************************** */
template<class Vector, class VectorSpace>
void
BDF1_TI<Vector, VectorSpace>::ReportStatistics_()
{
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::ostringstream oss;
    oss.flush();
    oss.setf(std::ios::scientific, std::ios::floatfield);

    oss << "TS:" << std::right << state_->seq;
    oss << " FS:" << state_->failed_solve;
    oss << " NS:" << state_->solve_itrs;
    oss << " PC:" << state_->pc_updates << " " << state_->pc_calls;
    oss << " LS:" << fn_->ReportStatistics();

    oss << " dt:";
    oss.precision(4);
    oss.width(10);
    oss << state_->hmin << " " << state_->hmax;

    *vo_->os() << oss.str() << std::endl;
  }
}

} // namespace Amanzi

#endif
