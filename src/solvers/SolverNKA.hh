/*
  This is the Nonlinear Solver component of the Amanzi code.

  Interface for using NKA as a solver.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_NKA_SOLVER_
#define AMANZI_NKA_SOLVER_

#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"

#include "SolverFnBase.hh"
#include "Solver.hh"
#include "NKA_Base.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector>
class SolverNKA : public Solver<Vector> {
 public:
  SolverNKA(Teuchos::ParameterList& plist,
            const Teuchos::RCP<SolverFnBase<Vector> >& fn,
            const Vector& initvec);

  virtual int Solve(const Teuchos::RCP<Vector>& u);

 protected:
  void Init_();

 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBase<Vector> > fn_;
  Teuchos::RCP<NKA_Base<Vector> > nka_;

  Teuchos::RCP<VerboseObject> vo_;

  double nka_tol_;
  int nka_dim_;

 private:
  double tol_, overflow_tol_;

  int max_itrs_, seq;
  int fun_calls_, pc_calls_;
  int pc_lag_, update_pc_calls_;
  int nka_lag_space_, nka_lag_iterations_;
  int max_error_growth_factor_, max_du_growth_factor_;
  int max_divergence_count_;
  int monitor_;
};


template<class Vector>
SolverNKA<Vector>::SolverNKA(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<SolverFnBase<Vector> >& fn,
                             const Vector& initvec) :
    plist_(plist), fn_(fn)
{
  Init_();

  // create the NKA space
  nka_ = Teuchos::rcp(new NKA_Base<Vector>(nka_dim_, nka_tol_, initvec));
}


template<class Vector>
void SolverNKA<Vector>::Init_()
{
  // NKA space control
  nka_dim_ = plist_.get<int>("max nka vectors", 10);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_ - 1);
  nka_tol_ = plist_.get<double>("nka vector tolerance", 0.05);

  tol_ = plist_.get<double>("nonlinear tolerance", 1.e-6);
  overflow_tol_ = plist_.get<double>("diverged tolerance", 1.0e10);
  max_itrs_ = plist_.get<int>("limit iterations", 100);
  max_du_growth_factor_ = plist_.get<double>("max du growth factor", 1.0e5);
  max_error_growth_factor_ = plist_.get<double>("max error growth factor", 1.0e5);
  max_divergence_count_ = plist_.get<int>("max divergent iterations", 3);
  // monitor_ = plist_.get<string>("convergence monitor", "monitor update");

  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject("AmanziSolvers::NKA", plist_));
}


template<class Vector>
int SolverNKA<Vector>::Solve(const Teuchos::RCP<Vector>& u) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // restart the nonlinear solver (flush its history)
  nka_->Restart();

  // initialize the iteration counter
  int itr(0);

  // create storage
  Teuchos::RCP<Vector> res = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du_tmp = Teuchos::rcp(new Vector(*u));
  du_tmp->PutScalar(0.0);

  // variables to monitor the progress of the nonlinear solver
  double error(0.0), previous_error(0.0), l2_error(0.);
  double l2_error_initial(0.);
  double du_norm(0.0), previous_du_norm(0.0);
  int divergence_count(0);

  // nonlinear solver main loop
  do {
    seq++;

    // Check for too many nonlinear iterations.
    if (itr > max_itrs_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Solve failed " << itr << " iterations (max), error = " << error << std::endl;

      return -1;
    }

    // update the preconditioner if necessary
    int errc(0);
    if (itr%(pc_lag_ + 1)==0) {
      update_pc_calls_++;
      fn_->UpdatePreconditioner(u);
    }

    // increment iteration counter
    itr++;

    // Evaluate the nonlinear function.
    fun_calls_++;
    fn_->Residual(u, res);

    // If monitoring the residual, check for convergence.
    if (monitor_ == SOLVER_MONITOR_RESIDUAL) {
      previous_error = error;
      error = fn_->ErrorNorm(u, res);
      res->Norm2(&l2_error);

      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << itr << ": error(res) = " << error << std::endl
                  << itr << ": L2 error(res) = " << l2_error << std::endl;

      // attempt to catch non-convergence early
      if (itr == 1) {
        l2_error_initial = l2_error;
      } else if (itr > 8) {
        if (l2_error > l2_error_initial) {
          if (vo_->os_OK(Teuchos::VERB_LOW))
            *vo_->os() << "Solve not converging, L2 error " << l2_error
                       << " > " << l2_error_initial << " (initial L2 error)"
                       << std::endl;
          return -1;
        }
      }

      if (error < tol_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solve succeeded: " << itr << " iterations, error = "
                     << error << std::endl;
        return itr;

      } else if (error > overflow_tol_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solve failed, error " << error << " > "
                     << overflow_tol_ << " (dtol)" << std::endl;
        return -1;

      } else if ((itr > 1) && (error > max_error_growth_factor_ * previous_error)) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solver is threatening to overflow, FAIL." << std::endl
                     << "  || error || = " << error << ", || error_prev || = "
                     << previous_error << std::endl;
        return -1;
      }
    }

    // Apply the preconditioner to the nonlinear residual.
    pc_calls_++;
    fn_->ApplyPreconditioner(res, du_tmp);

    // Calculate the accelerated correction.
    if (itr <= nka_lag_space_) {
      // lag the NKA space, just use the PC'd update
      *du = *du_tmp;
    } else {
      if (itr <= nka_lag_iterations_) {
        // lag NKA's iteration, but update the space with this Jacobian info
        nka_->Correction(*du_tmp, *du, du);
        *du = *du_tmp;
      } else {
        // take the standard NKA correction
        nka_->Correction(*du_tmp, *du, du);
      }
    }

    // Hack the correction
    bool hacked = fn_->ModifyCorrection(res, u, du);
    if (hacked) {
      // if we had to hack things, it't not unlikely that the Jacobian
      // information we have is crap.  Take the hacked correction, and restart
      // NKA to start building a new Jacobian space.
      nka_->Restart();
    }

    // Make sure that we do not diverge and cause numerical overflow.
    previous_du_norm = du_norm;
    du->NormInf(&du_norm);

    if ((itr > 1) && (du_norm > max_du_growth_factor_ * previous_du_norm)) {
      // try to recover by restarting NKA
      nka_->Restart();

      // ... this is the first invocation of nka_correction with an empty
      // nka space, so we call it withoug du_last, since there isn't one
      nka_->Correction(*du_tmp, *du);

      // re-check du
      du->NormInf(&du_norm);
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Solver is threatening to overflow, attempting to restart NKA."
                   << std::endl
                   << "  ||du|| = " << du_norm << ", ||du_prev||=" << previous_du_norm << std::endl;

      // If it fails again, give up.
      if ((itr > 1) && (du_norm > max_du_growth_factor_ * previous_du_norm)) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
           *vo_->os() << "Solver is threatening to overflow, FAIL." << std::endl
                      << "  ||du||=" << du_norm << ", ||du_prev||=" << previous_du_norm << std::endl;
        return -1;
      }
    }

    // Keep track of diverging iterations
    if (itr > 1 && du_norm >= previous_du_norm) {
      // The solver is threatening to diverge.
      ++divergence_count;

      // If it does not recover quickly, abort.
      if (divergence_count == max_divergence_count_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solver is diverging repeatedly, FAIL." << std::endl;
        return -1;
      }
    } else {
      divergence_count = 0;
    }

    // Next solution iterate and error estimate: u  = u - du
    u->Update(-1.0, *du, 1.0);
    fn_->ChangedSolution();

    // Monitor the PC'd residual
    if (monitor_ == SOLVER_MONITOR_PCED_RESIDUAL) {
      previous_error = error;
      error = fn_->ErrorNorm(u, du_tmp);
      du_tmp->Norm2(&l2_error);
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << itr << ": error(PC(res)) = " << error << std::endl
                   << itr << ": L2 error(PC(res)) = " << l2_error << std::endl;

      // Check for convergence.
      if (error < tol_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solve succeeded: " << itr << " iterations, error = "
                     << error << std::endl;
        return itr;
      } else if (error > overflow_tol_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solve failed, error " << error << " > "
                     << overflow_tol_ << " (dtol)" << std::endl;
        return -1;

      } else if ((itr > 1) && (error > max_error_growth_factor_ * previous_error)) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solver is threatening to overflow, FAIL." << std::endl
                     << "  ||error|| = " << error << ", ||error_prev||="
                     << previous_error << std::endl;
        return -1;
      }
    }

    // Monitor the NKA'd PC'd residual
    if (monitor_ == SOLVER_MONITOR_UPDATE) {
      previous_error = error;
      error = fn_->ErrorNorm(u, du);
      du->Norm2(&l2_error);

      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << itr << ": error=" << error
                   << ",  L2-error=" << l2_error << std::endl;

      // Check for convergence.
      if (error < tol_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solve succeeded: " << itr << " iterations, error = "
                     << error << std::endl;
        return itr;
      } else if (error > overflow_tol_) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solve failed, error " << error << " > "
                     << overflow_tol_ << " (dtol)" << std::endl;
        return -1;
      } else if ((itr > 1) && (error > max_error_growth_factor_ * previous_error)) {
        if (vo_->os_OK(Teuchos::VERB_LOW))
          *vo_->os() << "Solver is threatening to overflow, FAIL." << std::endl
                     << "  ||error||=" << error << ", ||error_prev||="
                     << previous_error << std::endl;
        return -1;
      }
    }
  } while (true);
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
