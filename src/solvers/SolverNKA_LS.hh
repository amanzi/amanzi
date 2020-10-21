/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! NKA nonlinear solver with a line-search based on a Brendt minimization algorithm.

/*!

Does NKA, then checks if that correction has reduced the residual by at least a
tolerance.  If not, this uses a Brendt minimization algorithm to try and find
an :math:`\alpha` that minimizes the reduction in the residual.

Note, this always monitors the residual.

.. _solver-nka-ls-spec:
.. admonition solver-nka-ls-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** Defines the required error
      tolerance. The error is calculated by a PK.

    * `"limit iterations`" ``[int]`` **20** Defines the maximum allowed number
      of iterations.

    * `"backtrack monitor`" ``[string]`` **monitor either** What norm is
      checked to determine whether backtracking has improved the residual or
      not?  One of `"monitor enorm`", `"monitor L2 residual`", or `'monitor
      either`"

    * `"backtrack tolerance`" ``[double]`` **0.** If the default update reduces
      the residual by at least this much, line search is not performed.

    * `"accuracy of line search minimum [bits]`" ``[int]`` **10** Convergence criteria on Brendt algorithm.

    * `"min valid alpha`" ``[double]`` **0** Lower bound on Brendt algorithm.

    * `"max valid alpha`" ``[double]`` **10.** Upper bound on Brendt algorithm.

    * `"max line search iterations`" ``[int]`` **10** Max iterations for the Brendt algorithm.

    * `"max nka vectors`" ``[int]`` **10** Defines the maximum number of
      consecutive vectors used for a local space.

    * `"nka vector tolerance`" ``[double]`` **0.05** Defines the minimum
      allowed orthogonality between vectors in the local space. If a new vector
      does not satisfy this requirement, the space is modified.

  
*/

#ifndef AMANZI_NKA_LINESEARCH_SOLVER_
#define AMANZI_NKA_LINESEARCH_SOLVER_

#include "boost/math/tools/minima.hpp"

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

template<class Vector, class VectorSpace>
class SolverNKA_LS : public Solver<Vector, VectorSpace> {
 public:
  SolverNKA_LS(Teuchos::ParameterList& plist) :
      plist_(plist) {};

  SolverNKA_LS(Teuchos::ParameterList& plist,
               const Teuchos::RCP<SolverFnBase<Vector> >& fn,
               const VectorSpace& map) :
      plist_(plist) {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
            const VectorSpace& map);

  int Solve(const Teuchos::RCP<Vector>& u) {
    returned_code_ = NKA_LS_(u);
    return (returned_code_ >= 0) ? 0 : 1;
  }

  // mutators
  void set_tolerance(double tol) { tol_ = tol; }
  void set_pc_lag(double pc_lag) { pc_lag_ = pc_lag; }
  virtual void set_db(const Teuchos::RCP<ResidualDebugger>& db) {
    db_ = db;
  }

  // access
  double tolerance() { return tol_; }
  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }
  int pc_calls() { return pc_calls_; }
  int pc_updates() { return pc_updates_; }
  int returned_code() { return returned_code_; }

 private:
  void Init_();
  int NKA_LS_(const Teuchos::RCP<Vector>& u);
  int NKA_ErrorControl_(double error, double previous_error, double l2_error);

 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBase<Vector> > fn_;
  Teuchos::RCP<NKA_Base<Vector, VectorSpace> > nka_;
  Teuchos::RCP<ResidualDebugger> db_;
  Teuchos::RCP<VerboseObject> vo_;

  // common parameters
  double tol_;
  int max_itrs_;
  int pc_lag_;

  // nka parameters
  double nka_tol_;
  int nka_dim_;

  // backtracking parameters
  double backtrack_atol_;
  double backtrack_rtol_;
  BacktrackMonitor bt_monitor_;
  
  // line search parameters
  int bits_;
  double min_alpha_;
  double max_alpha_;
  int max_ls_itrs_;

  // diagnostics
  double residual_;
  int num_itrs_;
  int returned_code_;
  int fun_calls_, pc_calls_, solve_calls_;
  int pc_updates_;


  // functor for minimization in boost
  struct Functor {
    Functor(const Teuchos::RCP<SolverFnBase<Vector> >& my_fn) :
        fn(my_fn) {}

    void setup(const Teuchos::RCP<Vector>& u_,
               const Teuchos::RCP<Vector>& u0_,
               const Teuchos::RCP<Vector>& du_) {
      u = u_;
      du = du_;
      u0 = u0_;
      if (r == Teuchos::null) {
        r = Teuchos::rcp(new Vector(*u));
      }
      *u0 = *u;
    }
    
    double operator()(double x) {
      *u = *u0;
      u->Update(-x, *du, 1.);
      fn->ChangedSolution();
      fn->Residual(u, r);
      return fn->ErrorNorm(u, r);
    }

    Teuchos::RCP<Vector> u,r,u0,du;
    Teuchos::RCP<SolverFnBase<Vector> > fn;
  };

};



/* ******************************************************************
* Public Init method.
****************************************************************** */
template<class Vector, class VectorSpace>
void
SolverNKA_LS<Vector,VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
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
template<class Vector, class VectorSpace>
void SolverNKA_LS<Vector, VectorSpace>::Init_()
{
  tol_ = plist_.get<double>("nonlinear tolerance", 1.e-6);
  max_itrs_ = plist_.get<int>("limit iterations", 20);
  residual_ = -1.0;

  // backtracking control
  std::string bt_monitor_string = plist_.get<std::string>("backtrack monitor",
          "monitor either");
  if (bt_monitor_string == "monitor enorm") {
    bt_monitor_ = BT_MONITOR_ENORM;
  } else if (bt_monitor_string == "monitor L2 residual") {
    bt_monitor_ = BT_MONITOR_L2;
  } else if (bt_monitor_string == "monitor either") {
    bt_monitor_ = BT_MONITOR_EITHER;
  } else {
    std::stringstream mstream;
    mstream << "Solver: Invalid backtrack monitor " << bt_monitor_string;
    Errors::Message m(mstream.str());
    Exceptions::amanzi_throw(m);
  }

  // backtracking control
  double backtrack_tol = plist_.get<double>("backtrack tolerance", 0.);
  backtrack_rtol_ = plist_.get<double>("backtrack relative tolerance", backtrack_tol);
  backtrack_atol_ = plist_.get<double>("backtrack absolute tolerance", backtrack_tol);

  // line search control
  bits_ = plist_.get<int>("accuracy of line search minimum [bits]", 10);
  min_alpha_ = plist_.get<double>("min valid alpha", 0.);
  max_alpha_ = plist_.get<double>("max valid alpha", 1.);
  max_ls_itrs_ = plist_.get<int>("max line search iterations", 10);

  // nka control
  nka_dim_ = plist_.get<int>("max nka vectors", 10);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_ - 1);
  nka_tol_ = plist_.get<double>("nka vector tolerance", 0.05);

  // diagnostics
  fun_calls_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;
  pc_lag_ = 0;
  solve_calls_ = 0;
  
  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject("Solver::NKA_LS", plist_));
}


/* ******************************************************************
* The body of NKA solver
****************************************************************** */
template<class Vector, class VectorSpace>
int SolverNKA_LS<Vector, VectorSpace>::NKA_LS_(const Teuchos::RCP<Vector>& u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  solve_calls_++;

  // restart the nonlinear solver (flush its history)
  nka_->Restart();

  // initialize the iteration and pc counters
  num_itrs_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;
  int db_write_iter = 0;

  // create storage
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du_tmp = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> u_prev = Teuchos::rcp(new Vector(*u));

  // variables to monitor the progress of the nonlinear solver
  double error(0.0), previous_error(0.0);
  double l2_error(0.0), previous_l2_error(0.0);

  // Evaluate the nonlinear function.
  fun_calls_++;
  fn_->Residual(u, r);

  // report error
  error = fn_->ErrorNorm(u, r);
  r->Norm2(&l2_error);
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << num_itrs_ << ": error(res) = " << error << std::endl
               << num_itrs_ << ": L2 error(res) = " << l2_error << std::endl;
  }

  // check if we have converged
  if (error < tol_) {
    residual_ = error;
    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os() << "Solve succeeded: " << num_itrs_ << " iterations, error = "
                 << error << std::endl;
    }
    return num_itrs_;
  }

  // set up the functor for minimization in line search
  Functor linesearch_func(fn_);
  linesearch_func.setup(u,u_prev,du);

  do {
    // Check for too many nonlinear iterations.
    if (num_itrs_ >= max_itrs_) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << "Solve reached maximum of iterations (" << num_itrs_ 
                   << ")  error=" << error << " terminating..." << std::endl;
      return SOLVER_MAX_ITERATIONS;
    }

    // Update the preconditioner if necessary.
    if (num_itrs_ % (pc_lag_ + 1) == 0) {
      pc_updates_++;
      fn_->UpdatePreconditioner(u);
    }

    // Increment iteration counter.
    num_itrs_++;

    // Apply the preconditioner to the nonlinear residual.
    pc_calls_++;
    fn_->ApplyPreconditioner(r, du_tmp);

    // Calculate the accelerated correction.
    nka_->Correction(*du_tmp, *du, du.ptr());

    // Hack the correction
    FnBaseDefs::ModifyCorrectionResult hacked = fn_->ModifyCorrection(r, u, du);
    if (hacked) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Hacked correction, restarting NKA." << std::endl;
      }
      // If we had to hack things, it't not unlikely that the Jacobian
      // information is crap. Take the hacked correction, and restart
      // NKA to start building a new Jacobian space.
      nka_->Restart();
    }
    
    // check the full correction
    bool good_iterate = false;
    bool admissible_iterate = false;
    *u_prev = *u;
    previous_error = error;
    previous_l2_error = l2_error;
    double alpha = 1.;
    u->Update(-alpha, *du, 1.);
    fn_->ChangedSolution();
    
    // Check the full correction
    if (fn_->IsAdmissible(u)) {
      admissible_iterate = true;

      // Evaluate the nonlinear function.
      fun_calls_++;
      fn_->Residual(u, r);

      // debugging
      db_->WriteVector<Vector>(db_write_iter++, *r, u.ptr(), du.ptr());
      
      // report error
      error = fn_->ErrorNorm(u, r);
      r->Norm2(&l2_error);
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << num_itrs_ << ": error(res) = " << error << std::endl
                   << num_itrs_ << ": L2 error(res) = " << l2_error << std::endl;
      }

      // check if we have converged
      if (error < tol_) {
        residual_ = error;
        if (vo_->os_OK(Teuchos::VERB_LOW)) {
          *vo_->os() << "Solve succeeded: " << num_itrs_ << " iterations, error = "
                     << error << std::endl;
        }
        return num_itrs_;
      }
      
      // Check if we have improved
      if (bt_monitor_ == BT_MONITOR_ENORM ||
          bt_monitor_ == BT_MONITOR_EITHER) {
        good_iterate |= error < previous_error * (1.+backtrack_rtol_) + backtrack_atol_;
      }
      if (bt_monitor_ == BT_MONITOR_L2 ||
          bt_monitor_ == BT_MONITOR_EITHER) {
        good_iterate |= l2_error < previous_l2_error * (1. + backtrack_rtol_) + backtrack_atol_;
      }
    }

    // If we haven't found a good correction yet, do a line search
    if ((!good_iterate) && (hacked != FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING)) {

      // find an admissible endpoint alpha, starting from ten times the full correction
      alpha = max_alpha_;
      *u = *u_prev;
      u->Update(-alpha, *du, 1.0);
      fn_->ChangedSolution();
      while (!fn_->IsAdmissible(u)) {
        alpha *= 0.3;
        *u = *u_prev;
        u->Update(-alpha, *du, 1.);
        fn_->ChangedSolution();
      }
      admissible_iterate = true;

      // minimize along the search path from min_alpha to endpoint
      double left = min_alpha_;
      boost::uintmax_t ls_itrs(max_ls_itrs_);
      std::pair<double,double> result = boost::math::tools::brent_find_minima(
          linesearch_func, left, alpha, bits_, ls_itrs);
      fun_calls_ += ls_itrs;
        
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "  Brent algorithm in: " << ls_itrs << " itrs (alpha="
                   << result.first << ") Error = " << result.second
                   << "(old error=" << previous_error << ")" << std::endl; 
      }
      alpha = result.first;
      error = result.second;

      // update the correction -- unclear if this is necessary -- TEST
      linesearch_func(result.first);
      fun_calls_++;

      // report error
      r->Norm2(&l2_error);
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << num_itrs_ << ": post-LS error(res) = " << error << std::endl
                   << num_itrs_ << ": post-LS L2 error(res) = " << l2_error << std::endl;
      }
      
      // check if we have converged
      if (error < tol_) {
        residual_ = error;
        if (vo_->os_OK(Teuchos::VERB_LOW)) {
          *vo_->os() << "Solve succeeded: " << num_itrs_ << " iterations, error = "
                     << error << std::endl;
        }
        return num_itrs_;
      }
      
      // Check if we have improved
      if (bt_monitor_ == BT_MONITOR_ENORM ||
          bt_monitor_ == BT_MONITOR_EITHER) {
        good_iterate |= error < previous_error * (1.+backtrack_rtol_) + backtrack_atol_;
      }
      if (bt_monitor_ == BT_MONITOR_L2 ||
          bt_monitor_ == BT_MONITOR_EITHER) {
        good_iterate |= l2_error < previous_l2_error * (1. + backtrack_rtol_) + backtrack_atol_;
      }

      // scale du to be the true correction
      du->Scale(alpha);

      // debugging
      db_->WriteVector<Vector>(db_write_iter++, *r, u.ptr(), du.ptr());
    }
    
    // Check for failure to find a good iteration
    if (!admissible_iterate) {
      // fail, bad search direction
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Solution iterate search direction not downhill, FAIL." << std::endl;
      }
      return SOLVER_BAD_SEARCH_DIRECTION;
    }

  } while (true);
}


}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
