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

template<class Vector, class VectorSpace>
class BDF1_TI {
 public:
  // Create the BDF Dae solver object, the nonlinear problem must be
  // defined in a class that derives from the virtual base class FnBase.
  BDF1_TI(BDFFnBase<Vector>& fn, Teuchos::ParameterList& plist,
          const Teuchos::RCP<const Vector>& initvector);

  // initializes the state
  void SetInitialState(const double h,
                       const Teuchos::RCP<Vector>& x,
                       const Teuchos::RCP<Vector>& xdot);

  // After a successful step, this method commits the new
  // solution to the solution history
  void CommitSolution(const double h, const Teuchos::RCP<Vector>& u);

  // computes a step
  bool TimeStep(double dt, double& dt_next, const Teuchos::RCP<Vector>& x);

  // Reset the memory of the time integrator
  void Reset();

  // returns the most recent time
  double time();

  // Report statistics
  void ReportStatistics(std::ostream&);

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

  AmanziSolvers::SolverFactory<Vector,VectorSpace> factory;
  solver_ = factory.Create(plist_);

  solver_->Init(solver_fn_, initvector->Map());
}


/* ******************************************************************
* Initialize miscaleneous parameters.
****************************************************************** */
template<class Vector,class VectorSpace>
void BDF1_TI<Vector,VectorSpace>::SetInitialState(const double t,
        const Teuchos::RCP<Vector>& x,
        const Teuchos::RCP<Vector>& xdot) {
  // set a clean initial state for when the time integrator is reinitialized
  state_->uhist->FlushHistory(t, *x, *xdot);
  state_->seq = 0;
  state_->pc_lag = 0;
}


/* ******************************************************************
* Record solution to the history.
****************************************************************** */
template<class Vector,class VectorSpace>
void BDF1_TI<Vector,VectorSpace>::CommitSolution(const double h, const Teuchos::RCP<Vector>& u) {
  double t = h + state_->uhist->MostRecentTime();

  // record the solution for later use when computing an initial guess
  // for the nonlinear solver
  state_->uhist->RecordSolution(t, *u);

  // record some information about this time step
  state_->hlast = h;
  state_->seq++;
  state_->freeze_count = std::max(0, state_->freeze_count-1);
  state_->hmin = std::min<double>(h, state_->hmin);
  state_->hmax = std::max<double>(h, state_->hmax);
}


/* ******************************************************************
* Returns most recent time.
****************************************************************** */
template<class Vector,class VectorSpace>
double BDF1_TI<Vector,VectorSpace>::time() {
  return state_->uhist->MostRecentTime();
}


/* ******************************************************************
* Implementation of implicit Euler time step.
****************************************************************** */
template<class Vector,class VectorSpace>
bool BDF1_TI<Vector,VectorSpace>::TimeStep(double dt, double& dt_next, const Teuchos::RCP<Vector>& u) {
  // initialize the output stream
  Teuchos::OSTab tab = vo_->getOSTab();

  // print some info about the time step
  double tlast = state_->uhist->MostRecentTime();
  double tnew = tlast + dt;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "step " << state_->seq << " T = " << tlast
               << " [sec]  dT = " << dt << std::endl;
    *vo_->os() << "preconditioner lag is " << state_->pc_lag
               << " out of " << state_->maxpclag << std::endl;
  }

  // u at the start of the time step
  Teuchos::RCP<Vector> u0 = Teuchos::rcp(new Vector(*u));

  // Predicted solution (initial value for the nonlinear solver)
  if (state_->extrapolate_guess) {
    if (state_->uhist->history_size() > 1) {
      state_->uhist->InterpolateSolution(tnew, *u);
      fn_->changed_solution();

      if (fn_->is_admissible(u)) {
	bool changed = fn_->ModifyPredictor(dt, u0, u);
	if (changed) fn_->changed_solution();
      } else {
	*u = *u0;
	fn_->changed_solution();
      }
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
      *vo_->os() << vo_->color("red") << "step failed with error code " << itr << vo_->reset() << std::endl;
      *u = *u0;
    }
  }

  // update the next timestep size
  dt_next = state_->ts_control->get_timestep(dt, itr);

  // update the preconditioner lag
  if (itr < 0) {
    state_->pc_lag = 0;
  } else {
    if (dt_next > dt) {
      state_->pc_lag = std::min(state_->pc_lag + 1, state_->maxpclag);
    } else if (dt_next < dt) {
      state_->pc_lag = std::max(state_->pc_lag - 1, 0);
    }
  }
  solver_->set_pc_lag(state_->pc_lag);

  // update performance statistics
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
    oss << std::left << state_.uhist->MostRecentTime();
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
void BDF1_TI<Vector,VectorSpace>::ReportStatistics(std::ostream& oss) {
  oss << "Report from BDF1 Time Integrator:" << std::endl;
  oss << "  total timesteps = " << state_->seq << std::endl;
  oss << "  total failed steps = " << state_->failed_bce << std::endl;
  oss << "  overall min dt = " << state_->hmin << std::endl;
  oss << "  overall max dt = " << state_->hmax << std::endl;
  oss << "--------------------------------------" << std::endl;
}

}  // namespace Amanzi

#endif
