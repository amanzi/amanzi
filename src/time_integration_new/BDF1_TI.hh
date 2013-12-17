#ifndef AMANZI_BDF1_TIME_INTEGRATOR_HH_
#define AMANZI_BDF1_TIME_INTEGRATOR_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "errors.hh"
#include "VerboseObject.hh"
#include "Solver.hh"
#include "SolverFactory.hh"

#include "BDF1_State.hh"
#include "SolutionHistory.hh"
#include "BDF1_SolverFnBase.hh"

namespace Amanzi {

template<class Vector,class VectorSpace>
class BDF1_TI {
 public:
  // Create the BDF Dae solver object, the nonlinear problem must be
  // defined in a class that derives from the virtual base class FnBase.
  BDF1_TI(BDFFnBase<Vector>& fn, Teuchos::ParameterList& plist,
          const Teuchos::RCP<const Vector>& initvector);

  // initializes the state
  void set_initial_state(const double h,
                         const Teuchos::RCP<Vector>& x,
                         const Teuchos::RCP<Vector>& xdot);

  // After a successful step, this method commits the new
  // solution to the solution history
  void commit_solution(const double h, const Teuchos::RCP<Vector>& u);

  // computes a step
  bool time_step(double dt, double& dt_next, const Teuchos::RCP<Vector>& x);

  // Reset the memory of the time integrator
  void reset();

  // returns the most recent time
  double time();

  // Report statistics
  void Report(std::ostream&);

 protected:
  void WriteSteppingStatistics_();

 protected:
  int mtries_;
  Teuchos::RCP<AmanziSolvers::Solver<Vector,VectorSpace> > solver_;
  Teuchos::RCP<BDF1_State<Vector> > state_;
  Teuchos::RCP<BDF1_SolverFnBase<Vector> > solver_fn_;
  Teuchos::RCP<BDFFnBase<Vector> > fn_;
  Teuchos::RCP<const Vector> initvector_;
  Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;
};


/* ******************************************************************
* Constructor
****************************************************************** */
template<class Vector,class VectorSpace>
BDF1_TI<Vector, VectorSpace>::BDF1_TI(BDFFnBase<Vector>& fn,
                     Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const Vector>& initvector) :
    plist_(plist), initvector_(initvector) {
  fn_ = Teuchos::rcpFromRef(fn);

  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject("TI::BDF1", plist_));

  // Create the state.
  state_ = Teuchos::rcp(new BDF1_State<Vector>());
  state_->InitializeFromPlist(plist_, initvector);

  // Set up the nonlinear solver
  // -- initialized the SolverFnBase interface
  solver_fn_ = Teuchos::rcp(new BDF1_SolverFnBase<Vector>(plist_, fn_));
  AmanziSolvers::SolverFactory<Vector,VectorSpace> fac;
  solver_ = fac.Create(plist_);
  solver_->Init(solver_fn_, initvector->Map());
}


/* ******************************************************************
* Initialize miscaleneous parameters.
****************************************************************** */
template<class Vector,class VectorSpace>
void BDF1_TI<Vector,VectorSpace>::set_initial_state(const double t,
        const Teuchos::RCP<Vector>& x,
        const Teuchos::RCP<Vector>& xdot) {
  // set a clean initial state for when the time integrator is reinitialized
  state_->uhist->flush_history(t, *x, *xdot);
  state_->seq = 0;
  state_->pc_lag = 0;
}


/* ******************************************************************
* Record solution to the history.
****************************************************************** */
template<class Vector,class VectorSpace>
void BDF1_TI<Vector,VectorSpace>::commit_solution(const double h, const Teuchos::RCP<Vector>& u) {
  double t = h + state_->uhist->most_recent_time();

  // record the solution for later use when computing an initial guess
  // for the nonlinear solver
  state_->uhist->record_solution(t, *u);

  // record some information about this time step
  state_->hlast = h;
  state_->seq++;
  state_->freeze_count = std::max<int>(0, state_->freeze_count-1);
  state_->hmin = std::min<double>(h, state_->hmin);
  state_->hmax = std::max<double>(h, state_->hmax);
}


/* ******************************************************************
* Returns most recent time.
****************************************************************** */
template<class Vector,class VectorSpace>
double BDF1_TI<Vector,VectorSpace>::time() {
  return state_->uhist->most_recent_time();
}


/* ******************************************************************
* Implementation of implicit Euler time step.
****************************************************************** */
template<class Vector,class VectorSpace>
bool BDF1_TI<Vector,VectorSpace>::time_step(double dt, double& dt_next, const Teuchos::RCP<Vector>& u) {
  // initialize the output stream
  Teuchos::OSTab tab = vo_->getOSTab();

  // print some info about the time step
  double tlast = state_->uhist->most_recent_time();
  double tnew = tlast + dt;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "step " << state_->seq + 1 << " T = " << tlast
               << " [sec]  dT = " << dt << std::endl;
  }

  // u at the start of the time step
  Teuchos::RCP<Vector> u0 = Teuchos::rcp(new Vector(*u));
  *u0 = *u;

  // Predicted solution (initial value for the nonlinear solver)
  if (state_->extrapolate_guess) {
    if ( state_->uhist->history_size() > 1) {
      state_->uhist->interpolate_solution(tnew,  *u);
      fn_->changed_solution();
      if (fn_->is_admissible(u)) {
	if (vo_->os_OK(Teuchos::VERB_HIGH))
	  *vo_->os() << "is admissible!" << std::endl;
	bool changed = fn_->modify_predictor(dt,u);
	if (changed) fn_->changed_solution();
      } else {
	*u = *u0;
	fn_->changed_solution();
      }
    } else {
      fn_->changed_solution();
      bool changed = fn_->modify_predictor(dt,u);
      if (changed) fn_->changed_solution();
    }
  }

  // Set up the solver fn
  solver_fn_->SetTimes(tlast, tnew);
  solver_fn_->SetPreviousTimeSolution(u0);

  // Solve the nonlinear BCE system.
  int itr;
  try {
    itr = solver_->Solve(u);
  } catch (const Errors::CutTimeStep& e) {
    itr = -1;
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    if (itr >= 0) {
      *vo_->os() << "success: " << solver_->num_itrs() << " nonlinear itrs" 
                 << " error=" << solver_->residual() << std::endl;
    } else {
      *vo_->os() << "failed with error code: " << itr << std::endl;
    }
  }

  // update the next timestep size
  dt_next = state_->ts_control->get_timestep(dt, itr);

  // update the pc lag
  if (itr < 0) {
    state_->pc_lag = 0;
  } else {
    if (dt_next > dt) {
      state_->pc_lag = std::min(state_->pc_lag + 1, state_->maxpclag);
    } else if (dt_next < dt) {
      state_->pc_lag = std::max(state_->pc_lag - 1, 0);
    }
  }

  // update performance monitors
  if (itr < 0) {
    state_->failed_bce++;
  } else {
    state_->hmax = std::max(state_->hmax, dt);
    state_->hmin = std::min(state_->hmin, dt);
  }
  return (itr < 0);
}


/* ******************************************************************
* Write statistics about the time step
****************************************************************** */
template<class Vector,class VectorSpace>
void BDF1_TI<Vector,VectorSpace>::WriteSteppingStatistics_() {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    std::ostringstream oss;

    oss.flush();
    oss.setf(std::ios::scientific, std::ios::floatfield);

    oss << "STEP=";
    oss.fill('0');
    oss.width(5);
    oss << std::right << state_.seq;
    oss << " T=";
    oss.precision(5);
    oss.width(11);
    oss << std::left << state_.uhist->most_recent_time();
    oss << " H=";
    oss.precision(5);
    oss.width(11);
    oss << std::left << state_.hlast;

    *vo_->os() << oss.str() << std::endl;
  }
}


/* ******************************************************************
* Report statistics.
****************************************************************** */
template<class Vector,class VectorSpace>
void BDF1_TI<Vector,VectorSpace>::Report(std::ostream& oss) {
  oss << "Report from BDF1 Time Integrator:" << std::endl;
  oss << "  total timesteps = " << state_->seq << std::endl;
  oss << "  total failed steps = " << state_->failed_bce << std::endl;
  oss << "  overall min dt = " << state_->hmin << std::endl;
  oss << "  overall max dt = " << state_->hmax << std::endl;
  oss << "--------------------------------------" << std::endl;
  //  solver_->Report(oss);
}

}  // namespace Amanzi

#endif
