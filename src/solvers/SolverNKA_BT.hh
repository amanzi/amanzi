/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for using NKA with backtracking as a solver.
*/

#ifndef AMANZI_NKA_BACKTRACKING_SOLVER_
#define AMANZI_NKA_BACKTRACKING_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"

#include "Solver.hh"
#include "SolverFnBase.hh"
#include "SolverDefs.hh"
#include "NKA_Base.hh"
#include "BackTracking.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverNKA_BT : public Solver<Vector, VectorSpace> {
 public:
  SolverNKA_BT(Teuchos::ParameterList& plist) : plist_(plist){};

  SolverNKA_BT(Teuchos::ParameterList& plist,
               const Teuchos::RCP<SolverFnBase<Vector>>& fn,
               const VectorSpace& map)
    : plist_(plist)
  {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn, const VectorSpace& map);

  int Solve(const Teuchos::RCP<Vector>& u)
  {
    returned_code_ = NKA_BT_(u);
    return (returned_code_ >= 0) ? 0 : 1;
  }

  // mutators
  void set_tolerance(double tol) { tol_ = tol; }
  void set_pc_lag(int pc_lag) { pc_lag_ = pc_lag; }

  // access
  double tolerance() { return tol_; }
  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }
  int pc_calls() { return pc_calls_; }
  int pc_updates() { return pc_updates_; }
  int returned_code() { return returned_code_; }

 private:
  void Init_();
  int NKA_BT_(const Teuchos::RCP<Vector>& u);
  int NKA_ErrorControl_(double error, double previous_error, double l2_error);

 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBase<Vector>> fn_;
  Teuchos::RCP<NKA_Base<Vector, VectorSpace>> nka_;

  Teuchos::RCP<VerboseObject> vo_;

  double nka_tol_;
  int nka_dim_;

 private:
  double tol_, overflow_tol_;

  int max_itrs_, num_itrs_, returned_code_;
  int fun_calls_, pc_calls_;
  int pc_lag_, pc_updates_;
  int nka_lag_space_, nka_lag_iterations_;
  int max_error_growth_factor_, max_du_growth_factor_;
  int max_divergence_count_;

  int backtrack_lag_, max_bt_;

  bool modify_correction_;
  double residual_; // defined by convergence criterion
  ConvergenceMonitor monitor_;
  int norm_type_;
};


/* ******************************************************************
* Public Init method.
****************************************************************** */
template <class Vector, class VectorSpace>
void
SolverNKA_BT<Vector, VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn,
                                        const VectorSpace& map)
{
  fn_ = fn;
  Init_();

  // Allocate the NKA space
  nka_ = Teuchos::rcp(new NKA_Base<Vector, VectorSpace>(nka_dim_, nka_tol_, map));
  nka_->Init(plist_);
}


/* ******************************************************************
* Initialization of the NKA solver.
****************************************************************** */
template <class Vector, class VectorSpace>
void
SolverNKA_BT<Vector, VectorSpace>::Init_()
{
  tol_ = plist_.get<double>("nonlinear tolerance", 1.e-6);
  overflow_tol_ = plist_.get<double>("diverged tolerance", 1.0e10);
  max_itrs_ = plist_.get<int>("limit iterations", 20);
  max_du_growth_factor_ = plist_.get<double>("max du growth factor", 1.0e5);
  max_error_growth_factor_ = plist_.get<double>("max error growth factor", 1.0e5);
  max_divergence_count_ = plist_.get<int>("max divergent iterations", 3);
  nka_lag_iterations_ = plist_.get<int>("lag iterations", 0);
  modify_correction_ = plist_.get<bool>("modify correction", false);

  std::string monitor_name = plist_.get<std::string>("monitor", "monitor update");
  ParseConvergenceCriteria(monitor_name, &monitor_, &norm_type_);

  nka_dim_ = plist_.get<int>("max nka vectors", 10);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_ - 1);
  nka_tol_ = plist_.get<double>("nka vector tolerance", 0.05);

  fun_calls_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;
  pc_lag_ = 0;
  nka_lag_space_ = 0;
  backtrack_lag_ = 0;

  residual_ = -1.0;

  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject("Solver::NKA_BT", plist_));
}


/* ******************************************************************
* The body of NKA solver
****************************************************************** */
template <class Vector, class VectorSpace>
int
SolverNKA_BT<Vector, VectorSpace>::NKA_BT_(const Teuchos::RCP<Vector>& u)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // restart the nonlinear solver (flush its history)
  nka_->Restart();

  // initialize the iteration and pc counters
  num_itrs_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;

  // create storage
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du_tmp = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> u_bt = Teuchos::rcp(new Vector(*u));

  // variables to monitor the progress of the nonlinear solver
  double error(0.0), previous_error(0.0), l2_error(0.0);
  double l2_error_initial(0.0);
  double du_norm(0.0), previous_du_norm(0.0);
  int divergence_count(0);

  do {
    // Check for too many nonlinear iterations.
    if (num_itrs_ >= max_itrs_) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << "Solve reached maximum of iterations (" << num_itrs_ << ")  error=" << error
                   << " terminating..." << std::endl;
      return SOLVER_MAX_ITERATIONS;
    }

    // Update the preconditioner if necessary.
    if (num_itrs_ % (pc_lag_ + 1) == 0) {
      pc_updates_++;
      fn_->UpdatePreconditioner(u);
    }

    // Increment iteration counter.
    num_itrs_++;

    // Evaluate the nonlinear function.
    fun_calls_++;
    fn_->Residual(u, r);

    // If monitoring the residual, check for convergence.
    if (monitor_ == SOLVER_MONITOR_RESIDUAL) {
      previous_error = error;
      error = fn_->ErrorNorm(u, r);
      residual_ = error;
      r->Norm2(&l2_error);

      // We attempt to catch non-convergence early.
      if (num_itrs_ == 1) {
        l2_error_initial = l2_error;
      } else if (num_itrs_ > 8) {
        if (l2_error > l2_error_initial) {
          if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
            *vo_->os() << "Solver stagnating, L2-error=" << l2_error << " > " << l2_error_initial
                       << " (initial L2-error)" << std::endl;
          return SOLVER_STAGNATING;
        }
      }

      int ierr = NKA_ErrorControl_(error, previous_error, l2_error);
      if (ierr == SOLVER_CONVERGED) return num_itrs_;
      if (ierr != SOLVER_CONTINUE) return ierr;
    }

    // Apply the preconditioner to the nonlinear residual.
    pc_calls_++;
    fn_->ApplyPreconditioner(r, du_tmp);

    // Calculate the accelerated correction.
    if (num_itrs_ <= nka_lag_space_) {
      // Lag the NKA space, just use the PC'd update.
      *du = *du_tmp;
    } else {
      if (num_itrs_ <= nka_lag_iterations_) {
        // Lag NKA's iteration, but update the space with this Jacobian info.
        nka_->Correction(*du_tmp, *du, du.ptr());
        *du = *du_tmp;
      } else {
        // Take the standard NKA correction.
        nka_->Correction(*du_tmp, *du, du.ptr());
      }
    }

    // Hack the correction
    if (modify_correction_) {
      bool hacked = fn_->ModifyCorrection(r, u, du);
      if (hacked) {
        // If we had to hack things, it't not unlikely that the Jacobian
        // information is crap. Take the hacked correction, and restart
        // NKA to start building a new Jacobian space.
        nka_->Restart();
      }
    }

    // Make sure that we do not diverge and cause numerical overflow.
    previous_du_norm = du_norm;
    du->NormInf(&du_norm);

    if ((num_itrs_ > 1) && (du_norm > max_du_growth_factor_ * previous_du_norm)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << "overflow: ||du||=" << du_norm << ", ||du_prev||=" << previous_du_norm
                   << std::endl
                   << "trying to restart NKA..." << std::endl;

      // Try to recover by restarting NKA.
      nka_->Restart();

      // This is the first invocation of nka_correction with an empty
      // nka space, so we call it withoug du_last, since there isn't one
      nka_->Correction(*du_tmp, *du);

      // Re-check du. If it fails again, give up.
      du->NormInf(&du_norm);

      if (du_norm > max_du_growth_factor_ * previous_du_norm) {
        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "terminating due to overflow ||du||=" << du_norm
                     << ", ||du_prev||=" << previous_du_norm << std::endl;
        return SOLVER_OVERFLOW;
      }
    }

    // Keep track of diverging iterations
    if (num_itrs_ > 1 && du_norm >= previous_du_norm) {
      divergence_count++;

      // If it does not recover quickly, abort.
      if (divergence_count == max_divergence_count_) {
        if (vo_->getVerbLevel() >= Teuchos::VERB_LOW)
          *vo_->os() << "Solver is diverging repeatedly, terminating..." << std::endl;
        return SOLVER_DIVERGING;
      }
    } else {
      divergence_count = 0;
    }

    // evaluate new step, backtrack if necessary
    if (num_itrs_ > backtrack_lag_) {
      previous_error = error;

      BackTracking<Vector> bt(fn_);
      bt.Bisection(u, du);

      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        *vo_->os() << bt.num_steps() << " backtracking steps,  ||r||: " << bt.initial_residual()
                   << " -> " << bt.final_residual() << std::endl;
      }
      nka_->Restart();
    }

    // Next solution iterate and error estimate: u  = u - du
    u->Update(-1.0, *du, 1.0);
    fn_->ChangedSolution();

    // Monitor the PC'd residual.
    if (monitor_ == SOLVER_MONITOR_PCED_RESIDUAL) {
      previous_error = error;
      error = fn_->ErrorNorm(u, du);
      residual_ = error;
      du_tmp->Norm2(&l2_error);

      int ierr = NKA_ErrorControl_(error, previous_error, l2_error);
      if (ierr == SOLVER_CONVERGED) return num_itrs_;
      if (ierr != SOLVER_CONTINUE) return ierr;
    }

    // Monitor the NKA'd PC'd residual.
    if (monitor_ == SOLVER_MONITOR_UPDATE) {
      previous_error = error;
      error = fn_->ErrorNorm(u, du);
      residual_ = error;
      du->Norm2(&l2_error);

      int ierr = NKA_ErrorControl_(error, previous_error, l2_error);
      if (ierr == SOLVER_CONVERGED) return num_itrs_;
      if (ierr != SOLVER_CONTINUE) return ierr;
    }
  } while (true);
}


/* ******************************************************************
* Internal convergence control.
****************************************************************** */
template <class Vector, class VectorSpace>
int
SolverNKA_BT<Vector, VectorSpace>::NKA_ErrorControl_(double error,
                                                     double previous_error,
                                                     double l2_error)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << num_itrs_ << ": error=" << error << "  L2-error=" << l2_error << std::endl;

  if (error < tol_) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Solver converged: " << num_itrs_ << " itrs, error=" << error << std::endl;
    return SOLVER_CONVERGED;
  } else if (overflow_tol_ > 0 && error > overflow_tol_) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
      *vo_->os() << "Solve failed, error " << error << " > " << overflow_tol_ << " (overflow)"
                 << std::endl;
    return SOLVER_OVERFLOW;
  } else if ((num_itrs_ > 1) && (error > max_error_growth_factor_ * previous_error)) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
      *vo_->os() << "Solver threatens to overflow, error " << error << " > " << previous_error
                 << " (previous error)" << std::endl;
    return SOLVER_OVERFLOW;
  }
  return SOLVER_CONTINUE;
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
