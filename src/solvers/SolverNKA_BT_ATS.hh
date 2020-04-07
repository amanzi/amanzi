/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for using NKA as a solver.
*/

#ifndef AMANZI_NKA_BT_ATS_SOLVER_
#define AMANZI_NKA_BT_ATS_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "ResidualDebugger.hh"

#include "errors.hh"
#include "FnBaseDefs.hh"
#include "Solver.hh"
#include "SolverFnBase.hh"
#include "SolverDefs.hh"
#include "NKA_Base.hh"
#include "AA_Base.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class SolverNKA_BT_ATS : public Solver<Vector, VectorSpace> {

 public:
  SolverNKA_BT_ATS(Teuchos::ParameterList& plist) :
      plist_(plist) {};

  SolverNKA_BT_ATS(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<SolverFnBase<Vector> >& fn,
                   const VectorSpace& map) :
      plist_(plist) {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
            const VectorSpace& map);

  virtual int Solve(const Teuchos::RCP<Vector>& u) {
    returned_code_ = NKA_BT_ATS_(u);
    return (returned_code_ >= 0) ? 0 : 1;
  }

  // mutators
  void set_tolerance(double tol) { tol_ = tol; }
  void set_pc_lag(double pc_lag) { pc_lag_ = pc_lag; }
  virtual void set_db(const Teuchos::RCP<ResidualDebugger>& db) {
    db_ = db;
  }

  // access
  double tolerance() {
    return tol_;
  }
  double residual() {
    return residual_;
  }
  int num_itrs() {
    return num_itrs_;
  }
  int pc_calls() {
    return pc_calls_;
  }
  int pc_updates() {
    return pc_updates_;
  }
  int returned_code() {
    return returned_code_;
  }

 private:
  void Init_();
  int NKA_BT_ATS_(const Teuchos::RCP<Vector>& u);
  int NKA_ErrorControl_(double error, double previous_error, double l2_error);

 private:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBase<Vector> > fn_;
  Teuchos::RCP<NKA_Base<Vector, VectorSpace> > nka_;
  Teuchos::RCP<AA_Base<Vector, VectorSpace> > aa_;
  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<ResidualDebugger> db_;

  double nka_tol_;
  int nka_dim_;

  double tol_, overflow_tol_;

  int max_itrs_, num_itrs_, returned_code_;
  int fun_calls_, pc_calls_, solve_calls_;
  int pc_lag_, pc_updates_;
  int nka_lag_iterations_;
  double nka_overflow_;

  bool fail_on_failed_backtrack_;
  double max_backtrack_;
  double max_total_backtrack_;
  double backtrack_lag_;
  double last_backtrack_iter_;
  double backtrack_factor_;
  double backtrack_atol_;
  double backtrack_rtol_;
  BacktrackMonitor bt_monitor_;

  double residual_;  // defined by convergence criterion

  bool use_aa_;
  bool use_nka_;
  double aa_beta_;
};


/* ******************************************************************
 * Public Init method.
 ****************************************************************** */
template<class Vector, class VectorSpace>
void
SolverNKA_BT_ATS<Vector,VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
        const VectorSpace& map)
{
  fn_ = fn;
  Init_();

  // Allocate the NKA space
  if (use_aa_) {
    aa_ = Teuchos::rcp(new AA_Base<Vector, VectorSpace>(nka_dim_, nka_tol_, aa_beta_, map));
    aa_->Init(plist_);
  } else {
    nka_ = Teuchos::rcp(new NKA_Base<Vector, VectorSpace>(nka_dim_, nka_tol_, map));
    nka_->Init(plist_);
  }

  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject(map.Comm(), "Solver::NKA_BT_ATS", plist_));
}


/* ******************************************************************
 * Initialization of the NKA solver.
 ****************************************************************** */
template<class Vector, class VectorSpace>
void SolverNKA_BT_ATS<Vector, VectorSpace>::Init_()
{
  tol_ = plist_.get<double>("nonlinear tolerance", 1.e-6);
  overflow_tol_ = plist_.get<double>("diverged tolerance", 1.0e10);
  max_itrs_ = plist_.get<int>("limit iterations", 20);

  nka_lag_iterations_ = plist_.get<int>("nka lag iterations", 0);
  nka_dim_ = plist_.get<int>("nka max vectors", 10);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_ - 1);
  nka_tol_ = plist_.get<double>("nka vector tolerance", 0.05);
  nka_overflow_ = plist_.get<double>("nka overflow protection", 1.e10);

  fun_calls_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;
  pc_lag_ = 0;
  solve_calls_ = 0;

  max_backtrack_ = plist_.get<int>("backtrack max steps", 10);
  max_total_backtrack_ = plist_.get<int>("backtrack max total steps", 1e6);
  backtrack_lag_ = plist_.get<int>("backtrack lag iterations", 0);
  last_backtrack_iter_ = plist_.get<int>("backtrack last iteration", 1e6);
  backtrack_factor_ = plist_.get<double>("backtrack factor", 0.5);

  double backtrack_tol = plist_.get<double>("backtrack tolerance", 0.);
  backtrack_rtol_ = plist_.get<double>("backtrack relative tolerance", backtrack_tol);
  backtrack_atol_ = plist_.get<double>("backtrack absolute tolerance", backtrack_tol);
  fail_on_failed_backtrack_ = plist_.get<bool>("backtrack fail on bad search direction", false);

  use_aa_ = plist_.get<bool>("Anderson mixing", false);
  if (use_aa_) {
    aa_beta_ = plist_.get<double>("relaxation parameter", 0.7);
  }

  std::string bt_monitor_string = plist_.get<std::string>("backtrack monitor",
          "monitor either");
  if (bt_monitor_string == "monitor enorm") {
    bt_monitor_ = BT_MONITOR_ENORM;
  } else if (bt_monitor_string == "monitor L2 residual") {
    bt_monitor_ = BT_MONITOR_L2;
  } else if (bt_monitor_string == "monitor either") {
    bt_monitor_ = BT_MONITOR_EITHER;
  } else {
    Errors::Message m;
    m << "Solver: Invalid backtrack monitor " << bt_monitor_string;
    Exceptions::amanzi_throw(m);
  }

  residual_ = -1.0;
}


/* ******************************************************************
 * The body of NKA solver
 ****************************************************************** */
template<class Vector, class VectorSpace>
int SolverNKA_BT_ATS<Vector, VectorSpace>::NKA_BT_ATS_(const Teuchos::RCP<Vector>& u) {
  solve_calls_++;

  // set the verbosity
  Teuchos::OSTab tab = vo_->getOSTab();

  // restart the nonlinear solver (flush its history)
  if (use_aa_) aa_->Restart();
  else nka_->Restart();

  // initialize the iteration and pc counters
  num_itrs_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;

  // create storage
  Teuchos::RCP<Vector> res = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du_nka = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du_pic = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> u_precorr = Teuchos::rcp(new Vector(*u));

  // variables to monitor the progress of the nonlinear solver
  double error(0.), previous_error(0.);
  double l2_error(0.), previous_l2_error(0.);
  bool nka_applied(false), nka_restarted(false);
  int nka_itr = 0;
  int total_backtrack = 0;
  int db_write_iter = 0;


  // Evaluate the nonlinear function.
  fun_calls_++;
  fn_->Residual(u, res);

  // Evaluate error
  error = fn_->ErrorNorm(u, res);
  db_->WriteVector<Vector>(db_write_iter++, *res, u.ptr(), du_nka.ptr());
  
  residual_ = error;
  res->Norm2(&l2_error);

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << num_itrs_ << ": error(res) = " << error << std::endl
               << num_itrs_ << ": L2 error(res) = " << l2_error << std::endl;
  }

  // Check divergence
  if (error > overflow_tol_) {
    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os() << "Solve failed, error " << error << " > " << overflow_tol_
                 << " (dtol)" << std::endl;
    }
    return SOLVER_DIVERGING;
  }

  // Check converged
  if (error < tol_) {
    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os() << "Solve succeeded: " << num_itrs_ << " iterations, error = "
                 << error << std::endl;
    }
    return num_itrs_;
  }

  // nonlinear solver main loop
  do {
    // increment iteration counter
    num_itrs_++;

    // evaluate precon at the beginning of the method
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Updating preconditioner." << std::endl;
    }
    pc_updates_++;
    fn_->UpdatePreconditioner(u);

    // Apply the preconditioner to the nonlinear residual.
    pc_calls_++;
    du_pic->PutScalar(0.);
    fn_->ApplyPreconditioner(res, du_pic);
    double du_pic_norm = 0;
    du_pic->NormInf(&du_pic_norm);

    if (nka_restarted) {
      // NKA was working, but failed.  Reset the iteration counter.
      nka_itr = 0;
      nka_restarted = false;
    }

    FnBaseDefs::ModifyCorrectionResult hacked = FnBaseDefs::CORRECTION_NOT_MODIFIED;
    if (num_itrs_ > nka_lag_iterations_ &&
        du_pic_norm < nka_overflow_) {
      // Calculate the accelerated correction.
      if (use_aa_) aa_->Correction(*du_pic, *du_nka, u.ptr());
      else nka_->Correction(*du_pic, *du_nka);
      nka_applied = true;
      nka_itr++;

      // Hack the correction
      hacked = fn_->ModifyCorrection(res, u, du_nka);
      if (hacked) {
        // if we had to hack things, it is likely that the Jacobian
        // information we have is, for the next iteration, not very accurate.
        // Take the hacked correction, and restart NKA to start building a new
        // Jacobian space.
        if (vo_->os_OK(Teuchos::VERB_HIGH)) {
          *vo_->os() << "Restarting NKA, correction modified." << std::endl;
        }
        if (use_aa_) aa_->Restart();
        else nka_->Restart();
        nka_restarted = true;
      }
    } else {
      // Do not calculate the accelarated correction.
      nka_applied = false;
      // Hack the Picard update
      hacked = fn_->ModifyCorrection(res, u, du_pic);

      if (du_pic_norm >= nka_overflow_) {
        if (vo_->os_OK(Teuchos::VERB_HIGH)) {
          *vo_->os() << "Restarting NKA, too large of a correction (overflow)." << std::endl;
        }
        if (use_aa_) aa_->Restart();
        else nka_->Restart();
        nka_restarted = true;
        
      }
    }

    // potentially backtrack
    bool admitted_iterate = false;
    if (num_itrs_ > backtrack_lag_ && num_itrs_ < last_backtrack_iter_ &&
        hacked != FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING) {
      bool good_step = false;
      *u_precorr = *u;
      previous_error = error;
      previous_l2_error = l2_error;

      if (nka_applied) {
        // check the NKA updated norm
        u->Update(-1, *du_nka, 1.);
        fn_->ChangedSolution();

        // Check admissibility of the iterate
        admitted_iterate = false;
        if (fn_->IsAdmissible(u)) {
          admitted_iterate = true;

          // Evaluate the nonlinear function.
          fun_calls_++;
          fn_->Residual(u, res);

          // Evalute error
          error = fn_->ErrorNorm(u, res);
          db_->WriteVector<Vector>(db_write_iter++, *res, u.ptr(), du_nka.ptr());
          
          residual_ = error;
          res->Norm2(&l2_error);
          if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
            *vo_->os() << num_itrs_ << ": NKA "
                       << ": error(res) = " << error << std::endl
                       << num_itrs_ << ": NKA "
                       << ": L2 error(res) = " << l2_error << std::endl;
          }

          // Check if we have improved
          if (bt_monitor_ == BT_MONITOR_ENORM ||
              bt_monitor_ == BT_MONITOR_EITHER) {
            good_step |= error < previous_error * (1.+backtrack_rtol_) + backtrack_atol_;
          }
          if (bt_monitor_ == BT_MONITOR_L2 ||
              bt_monitor_ == BT_MONITOR_EITHER) {
            good_step |= l2_error < previous_l2_error * (1. + backtrack_rtol_) + backtrack_atol_;
          }
        } // IsAdmissible()

        // if NKA did not improve the error, toss Jacobian info
        if (!good_step) {
          if (vo_->os_OK(Teuchos::VERB_HIGH)) {
            *vo_->os() << "Restarting NKA, NKA step does not improve error or resulted in inadmissible solution, on NKA itr = "
                       << nka_itr << std::endl;
          }
          if (use_aa_) aa_->Restart();
          else nka_->Restart();
          nka_restarted = true;

          // Tried NKA, but failed, so will now use Picard.  The Picard update
          // needs to be hacked (since it was not hacked in case NKA worked),
          // unless this is NKA itr 1, in which case the NKA update IS the
          // Picard update, so we can just copy over the NKA hacked version.
          if (nka_itr == 1) {
            *du_pic = *du_nka;
          } else {
            hacked = fn_->ModifyCorrection(res, u, du_pic);

            if (hacked == FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING) {
              // no backtracking, just use this correction, checking admissibility
              u->Update(-1., *du_pic, 1.);
              fn_->ChangedSolution();
              admitted_iterate = false;
              if (fn_->IsAdmissible(u)) {
                admitted_iterate = true;

                // Evaluate the nonlinear function.
                fun_calls_++;
                fn_->Residual(u, res);

                // Evalute error
                error = fn_->ErrorNorm(u, res);
                db_->WriteVector<Vector>(db_write_iter++, *res, u.ptr(), du_pic.ptr());

                residual_ = error;
                res->Norm2(&l2_error);
                if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
                  *vo_->os() << num_itrs_ << ": PIC "
                             << ": error(res) = " << error << std::endl
                             << num_itrs_ << ": PIC "
                             << ": L2 error(res) = " << l2_error << std::endl;
                }
              }
              good_step = true;
            }
          }
        }
      }

      if (!good_step) {
        // Backtrack on the Picard iterate
        bool done_backtracking = false;
        int n_backtrack = 0;
        double backtrack_alpha = 1.;

        if (nka_applied && nka_itr == 1) {
          // In nka_itr == 1, NKA is the Picard update.  We have already tried
          // it and it failed.  No need to try it a second time.
          n_backtrack++;
          backtrack_alpha *= backtrack_factor_;
          done_backtracking = n_backtrack > max_backtrack_;
        }

        while (!done_backtracking) {
          total_backtrack++;

          // apply the correction
          *u = *u_precorr;
          u->Update(-backtrack_alpha, *du_pic, 1.);
          fn_->ChangedSolution();

          // Check admissibility of the iterate
          admitted_iterate = false;
          if (fn_->IsAdmissible(u)) {
            admitted_iterate = true;

            // Evaluate the nonlinear function.
            fun_calls_++;
            fn_->Residual(u, res);

            // Evalute error
            error = fn_->ErrorNorm(u, res);
            db_->WriteVector<Vector>(db_write_iter++, *res, u.ptr(), du_pic.ptr());
            
            residual_ = error;
            res->Norm2(&l2_error);
            if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
              *vo_->os() << num_itrs_ << ": backtrack " << n_backtrack
                         << ": error(res) = " << error << std::endl
                         << num_itrs_ << ": backtrack " << n_backtrack
                         << ": L2 error(res) = " << l2_error << std::endl;
            }

            // Check if we have improved
            if (bt_monitor_ == BT_MONITOR_ENORM ||
                bt_monitor_ == BT_MONITOR_EITHER) {
              good_step |= error < previous_error * (1.+backtrack_rtol_) + backtrack_atol_;
            }
            if (bt_monitor_ == BT_MONITOR_L2 ||
                bt_monitor_ == BT_MONITOR_EITHER) {
              good_step |= l2_error < previous_l2_error * (1. + backtrack_rtol_) + backtrack_atol_;
            }
          } // IsAdmissible()

          if (!good_step) {
            // Not good, either because not admissible or not error-decreasing.
            n_backtrack++;
            backtrack_alpha *= backtrack_factor_;
          }

          // Check for done
          if (fail_on_failed_backtrack_ && n_backtrack > max_backtrack_) {
            // fail, bad search direction
            if (vo_->os_OK(Teuchos::VERB_LOW)) {
              *vo_->os() << "Solution iterate search direction not downhill, FAIL." << std::endl;
            }
            return SOLVER_BAD_SEARCH_DIRECTION;
          }
          done_backtracking = good_step || n_backtrack > max_backtrack_;
        } // backtrack loop
      } // backtracking

    } else {
      // not backtracking, either due to BT Lag or past the last BT iteration
      // apply the correction
      if (nka_applied) {
        u->Update(-1., *du_nka, 1.);
        fn_->ChangedSolution();
      } else {
        u->Update(-1., *du_pic, 1.);
        fn_->ChangedSolution();
      }

      // check correction admissibility
      admitted_iterate = false;
      if (fn_->IsAdmissible(u)) {
        admitted_iterate = true;

        // Evaluate the nonlinear function.
        fun_calls_++;
        fn_->Residual(u, res);

        // Evalute error
        error = fn_->ErrorNorm(u, res);
        if (nka_applied) {
          db_->WriteVector<Vector>(db_write_iter++, *res, u.ptr(), du_nka.ptr());
        } else {
          db_->WriteVector<Vector>(db_write_iter++, *res, u.ptr(), du_pic.ptr());
        }
        
        residual_ = error;
        res->Norm2(&l2_error);
        if (vo_->os_OK(Teuchos::VERB_LOW)) {
          *vo_->os() << num_itrs_ << (nka_applied ? ": NKA " : ": PIC ")
                     << ": error(res) = " << error << std::endl
                     << num_itrs_ << (nka_applied ? ": NKA " : ": PIC ")
                     << ": L2 error(res) = " << l2_error << std::endl;
        }
      }
    } // non-bt fork

    // Check final admissibility
    if (!admitted_iterate) {
      // final update was never evaluated, not admissible
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Solution iterate is not admissible, FAIL." << std::endl;
      }
      return SOLVER_INADMISSIBLE_SOLUTION;
    }

    // Check divergence
    if (error > overflow_tol_) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Solve failed, error " << error << " > " << overflow_tol_
                   << " (dtol)" << std::endl;
      }
      return SOLVER_DIVERGING;
    }

    // Check converged
    if (error < tol_) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Solve succeeded: " << num_itrs_ << " iterations, error = "
                   << error << std::endl;
      }
      return num_itrs_;
    }

    // Check for too many nonlinear iterations.
    if (num_itrs_ > max_itrs_) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Solve failed " << num_itrs_ << " iterations (max), error = "
                   << error << std::endl;
      }
      return SOLVER_MAX_ITERATIONS;
    }

    // Check for too many backtrack steps.
    if (total_backtrack > max_total_backtrack_) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Solve failed, backtracked " << total_backtrack << " times, giving up with error = "
                   << error << std::endl;
      }
      return SOLVER_STAGNATING;
    }
  } while (true);
}

} // namespace
} // namespace
#endif
