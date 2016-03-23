/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for using Backtracking line search as a solver.

  From PETSc SNES type BT, which in turn is from Numerical Methods for
  Unconstrained Optimization and Nonlinear Equations by Dennis & Schnabel, pg
  325.
*/

#ifndef AMANZI_BT_SOLVER_
#define AMANZI_BT_SOLVER_

#include "boost/math/tools/minima.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "ResidualDebugger.hh"

#include "Solver.hh"
#include "SolverFnBase.hh"
#include "SolverDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class SolverBT : public Solver<Vector,VectorSpace> {
 public:
  SolverBT(Teuchos::ParameterList& plist) :
      plist_(plist) {};

  SolverBT(Teuchos::ParameterList& plist,
           const Teuchos::RCP<SolverFnBase<Vector> >& fn,
           const VectorSpace& map) :
      plist_(plist) {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
            const VectorSpace& map);

  int Solve(const Teuchos::RCP<Vector>& u) {
    returned_code_ = BT_(u);
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
  int BT_(const Teuchos::RCP<Vector>& u);
  int BT_ErrorControl_(double error, double previous_error, double l2_error);
  

  struct Functor {
    Functor(const Teuchos::RCP<SolverFnBase<Vector> >& fn_) :
        fn(fn_) {}

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
      fn_calls = 0;
      residual_x = 0;
    }
    
    double operator()(double x) {
      fn_calls++;
      *u = *u0;
      u->Update(-x, *du, 1.);
      fn->Residual(u, r);
      residual_x = x;
      return fn->ErrorNorm(u, r);
    }

    Teuchos::RCP<Vector> u,r,u0,du;
    Teuchos::RCP<SolverFnBase<Vector> > fn;
    int fn_calls;
    double residual_x;
  };

  
 private:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBase<Vector> > fn_;

  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<ResidualDebugger> db_;

 private:
  double tol_, overflow_tol_, overflow_l2_tol_, overflow_pc_tol_, overflow_r_tol_;

  int max_itrs_, num_itrs_, returned_code_;
  int fun_calls_, pc_calls_;
  int pc_lag_, pc_updates_;
  int nka_lag_space_, nka_lag_iterations_;
  int max_error_growth_factor_, max_du_growth_factor_;
  int max_divergence_count_;

  int bits_;
  double alpha_;
  int max_ls_itrs_;
  
  bool modify_correction_;
  double residual_;  // defined by convergence criterion
  ConvergenceMonitor monitor_;
};



/* ******************************************************************
* Public Init method.
****************************************************************** */
template<class Vector, class VectorSpace>
void
SolverBT<Vector,VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
        const VectorSpace& map)
{
  fn_ = fn;
  Init_();
}


/* ******************************************************************
* Initialization of the NKA solver.
****************************************************************** */
template<class Vector, class VectorSpace>
void SolverBT<Vector, VectorSpace>::Init_()
{
  tol_ = plist_.get<double>("nonlinear tolerance", 1.e-6);
  overflow_tol_ = plist_.get<double>("diverged tolerance", 1.0e10);
  overflow_l2_tol_ = plist_.get<double>("diverged l2 tolerance", 1.0e10);
  overflow_pc_tol_ = plist_.get<double>("diverged pc tolerance", 1.0e10);
  overflow_r_tol_ = plist_.get<double>("diverged residual tolerance", 1.0e10);
  max_itrs_ = plist_.get<int>("limit iterations", 20);
  max_du_growth_factor_ = plist_.get<double>("max du growth factor", 1.0e5);
  max_error_growth_factor_ = plist_.get<double>("max error growth factor", 1.0e5);
  max_divergence_count_ = plist_.get<int>("max divergent iterations", 3);
  modify_correction_ = plist_.get<bool>("modify correction", false);

  bits_ = plist_.get<int>("minimization bits", 10);
  alpha_ = plist_.get<double>("minimum reduction factor", 0.00001);
  max_ls_itrs_ = plist_.get<int>("max line search iterations", 10);

  fun_calls_ = 0;
  pc_calls_ = 0;
  pc_updates_ = 0;
  pc_lag_ = 0;
  nka_lag_space_ = 0;

  residual_ = -1.0;

  
  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject("Solver::BT", plist_));
}


template<class Vector, class VectorSpace>
int
SolverBT<Vector,VectorSpace>::BT_(const Teuchos::RCP<Vector>& u)
{
  
  // create storage
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> r_end = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> du = Teuchos::rcp(new Vector(*u));
  Teuchos::RCP<Vector> u0 = Teuchos::rcp(new Vector(*u));

  // variables to monitor the progress of the nonlinear solver
  double error(0.0), previous_error(0.0);
  double l2_error(0.0), l2_error_initial(0.0);
  double du_norm(0.0), du_norm_initial(0.), previous_du_norm(0.0), r_norm_initial;
  int divergence_count(0);
  int prec_error;
  int db_write_iter = 0;

  num_itrs_ = 0;

  // Evaluate the nonlinear function.
  fun_calls_++;
  fn_->Residual(u, r);
  db_->WriteVector<Vector>(db_write_iter++, *r, u.ptr(), du.ptr());

  // If monitoring the residual, check for convergence.
  error = fn_->ErrorNorm(u, r);
  previous_error = error;
  residual_ = error;
  r->Norm2(&l2_error);
  l2_error_initial = l2_error;

  int ierr = BT_ErrorControl_(error, previous_error, l2_error);
  if (ierr == SOLVER_CONVERGED) return num_itrs_;
  if (ierr != SOLVER_CONTINUE) return ierr;

  // set up the functor for minimization in the line search
  Functor linesearch_func(fn_);
  linesearch_func.setup(u,u0,du);

  // loop til convergence or failure
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
    prec_error = fn_->ApplyPreconditioner(r, du);

    // Hack the correction
    if (modify_correction_) fn_->ModifyCorrection(r, u, du);

    // find an admissible endpoint, starting from the full correction
    double endpoint = 1.;
    *u0 = *u;
    u->Update(-endpoint, *du, 1.0);
    while (!fn_->IsAdmissible(u)) {
      endpoint *= 0.1;
      *u = *u0;
      u->Update(-endpoint, *du, 1.);
    }

    // minimize
    linesearch_func.fn_calls = 0;
    double left = 0.0;
    boost::uintmax_t max_itrs(max_ls_itrs_);
    std::pair<double,double> result = boost::math::tools::brent_find_minima(
        linesearch_func, left, endpoint, bits_, max_itrs);
    fun_calls_ += linesearch_func.fn_calls;

    std::cout << "   Brent found:" << result.first << ", " << result.second << std::endl;

    // check for a minimization value at the start
    if (std::abs(result.second - error) < std::max(1.e-10, error*alpha_)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << "Searching in this direction resulted in error=" << error << " which is not a sufficient reduction, indicating a bad search direction..." << std::endl;
      return SOLVER_BAD_SEARCH_DIRECTION;
    }

    // update the correction
    *u = *u0;
    u->Update(-result.first, *du, 1.);

    if (linesearch_func.residual_x != result.first) {
      // get the residual again
      fn_->Residual(u, r);
    }

    // test convergence
    previous_error = error;
    error = result.second;
    residual_ = error;
    r->Norm2(&l2_error);

    int ierr = BT_ErrorControl_(error, previous_error, l2_error);
    if (ierr == SOLVER_CONVERGED) return num_itrs_;
    if (ierr != SOLVER_CONTINUE) return ierr;
  } while(true);
}    


/* ******************************************************************
* Internal convergence control.
****************************************************************** */
template<class Vector, class VectorSpace>
int SolverBT<Vector, VectorSpace>::BT_ErrorControl_(
   double error, double previous_error, double l2_error)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) 
    *vo_->os() << num_itrs_ << ": error=" << error << "  L2-error=" << l2_error << std::endl;

  if (error < tol_) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) 
      *vo_->os() << "Solver converged: " << num_itrs_ << " itrs, error=" << error << std::endl;
    return SOLVER_CONVERGED;
  } else if (error > overflow_tol_) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) 
      *vo_->os() << "Solve failed, error " << error << " > "
                 << overflow_tol_ << " (overflow)" << std::endl;
    return SOLVER_OVERFLOW;
  } else if ((num_itrs_ > 1) && (error > max_error_growth_factor_ * previous_error)) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) 
      *vo_->os() << "Solver threatens to overflow, error " << error << " > "
                 << previous_error << " (previous error)" << std::endl;
    return SOLVER_OVERFLOW;
  }
  return SOLVER_CONTINUE;
}


} // namespace Amanzi 
} // namespace AmanziSolvers


#endif
