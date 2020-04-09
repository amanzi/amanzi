/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! A framework for nonlinear solvers that use the same basic structure.

/*!

The majority of nonlinear solvers use a common structure to implement their
solve process.  This reduces that to a few common callbacks, making it easier
to share common implementation here.

Implementations must implement the CalculateUpdate_() method, and may
re-implement the Restart_() and Init() methods.

 * `"monitor`" ``[string]`` **monitor l2 update** Defines what is meant by "error."  One of:
        - `"monitor l2 update`" The L2 norm of the correction, du.
        - `"monitor linf update`" The Linf norm of the correction, du
        - `"monitor update`" The functional's ErrorNorm() of the correction, du.
        - `"monitor l2 residual`" The L2 norm of the functional's residual.
        - `"monitor linf residual`" The Linf norm of the functional's residual.
        - `"monitor residual`" The functional's ErrorNorm() of the functional's residual.
      Hereafter this quantity is called the error; all tolerances, etc, are with respect to this.

 * `"nonlinear tolerance`" ``[double]`` **1.e-6** Solve is successful when the
       error is less than this tolerance.

 * `"limit iterations`" ``[int]`` **50** Max number of iterations before solver fails.

 * `"max divergent iterations`" ``[int]`` **3** Max number of iterations where
        the error increases before failure.

 * `"diverged tolerance`" ``[double]`` **1.e10** Max error allowed without failure.

 * `"max error growth factor`" ``[double]`` **1.e5** Max ratio of error to
        previous iteration's error allowed without failure.

 * `"max iterations before stagnation`" ``[int]`` **8** Number of iterations
        after which the error must have decreased from the initial error
        without failing.
 
 * `"modify correction`" ``[bool]`` **false** Call the PK's ModifyCorrection() method.

 * `"preconditioner lag iterations`" ``[int]`` **0** Number of iterations to
        lag before updating the preconditioner.

*/

#ifndef AMANZI_SOLVER_DEFAULT_
#define AMANZI_SOLVER_DEFAULT_

#include "Teuchos_RCP.hpp"

#include "ResidualDebugger.hh"
#include "Solver.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverDefault : public Solver<Vector,VectorSpace> {
 public:
  SolverDefault(Teuchos::ParameterList& plist);
  
  virtual void Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn,
                    const Teuchos::RCP<const VectorSpace>& map) override
  { fn_ = fn; }

  // Returns 0 if success, 1 if failure.
  virtual int Solve(const Teuchos::RCP<Vector>& u) override;
  
  // mutators
  virtual void set_tolerance(double tol) override { tolerance_ = tol; }
  virtual void set_pc_lag(int pc_lag) override { pc_lag_ = pc_lag; }
  virtual void set_db(const Teuchos::RCP<ResidualDebugger>& db) override
  { db_ = db; }

  // accessors
  // Note, what the error is is dependent upon the MonitorType and the
  // MonitorNorm
  virtual double error() const override { return error_; }

  // Tolerance to compare to error.
  virtual double tolerance() const override { return tolerance_; }

  // The L2 norm of the residual
  virtual double residual() const override { return residual_; }

  // Number of nonlinear iterations
  virtual int num_iterations() const override { return num_itrs_; }

  // See SolverDefs.hh MonitorStatus definition, but positive values indicate
  // number of iterations convergence was achieved in while negative numbers
  // indicate an error.
  virtual int returned_code() const override {
    if (status_ == MonitorStatus::CONVERGED) return num_iterations();
    else return (int) status_;
  }

  // Number of preconditioner ApplyInverse calls (this solve)
  virtual int pc_calls() const override { return pc_calls_; }

  // Number of times predonditioner was updated (this solve)
  virtual int pc_updates() const override { return pc_updates_; }

  // Number of times nonlinear residual function was called (this solve)
  virtual int function_calls() const override { return function_calls_; }


 protected:
  virtual void Restart_() {}

  // Do any modifications specific to the method.
  //
  // This may be, for instance, apply an accelerator, doing line search, etc,
  // or it may simply call the PK's ModifyCorrection() (it probably, at the
  // least, should do this).
  //
  // NOTE: It is up to the method to ensure that, after this call, if status is
  // CONTINUE, that du is valid, not diverging, and good to directly apply.
  //
  // NOTE: du is taken by non-const reference -- the pointer MAY be reseated.
  virtual std::pair<MonitorStatus,double> ModifyCorrection_(Teuchos::RCP<Vector>& r,
          const Teuchos::RCP<Vector>& u, Teuchos::RCP<Vector>& du) = 0;
  virtual MonitorStatus MonitorError_(double error, double previous_error,
          double l2_error, double l2_error_initial, int& divergence_count);
  
 protected:
  // public interface
  double tolerance_;
  double residual_;
  double error_;
  int num_itrs_;
  int returned_code_;
  int pc_calls_;
  int pc_updates_;
  int function_calls_;

  // functional 
  Teuchos::RCP<SolverFnBase<Vector>> fn_;

  // convergence control
  Monitor monitor_type_;
  MonitorStatus status_;
  MonitorNorm norm_type_;
  bool modify_correction_;
  int pc_lag_;

  // reasons to crash
  double diverged_tol_;
  double max_error_growth_factor_;
  int max_divergence_count_;
  int max_itrs_stagnation_;
  int max_itrs_;

  // debugging and verbose object
  Teuchos::RCP<ResidualDebugger> db_;
  Teuchos::RCP<VerboseObject> vo_;
};


// Implementation
template <class Vector, class VectorSpace>
SolverDefault<Vector,VectorSpace>::SolverDefault(Teuchos::ParameterList& plist)
    : function_calls_(0),
      pc_calls_(0),
      pc_updates_(0),
      pc_lag_(0),
      residual_(-1.0)
{
  // control parameters
  modify_correction_ = plist.get<bool>("modify correction", false);
  pc_lag_ = plist.get<int>("preconditioner lag iterations", 0);

  // converged/diverged control parameters
  tolerance_ = plist.get<double>("nonlinear tolerance", 1.e-6);
  max_itrs_ = plist.get<int>("limit iterations", 50);
  max_divergence_count_ = plist.get<int>("max divergent iterations", 3);
  diverged_tol_ = plist.get<double>("diverged tolerance", 1.0e10);
  max_error_growth_factor_ = plist.get<double>("max error growth factor", 1.0e5);
  max_itrs_stagnation_ = plist.get<int>("max iterations before stagnation", 8);

  // what do we monitor and how do we monitor it?
  std::string monitor_name = plist.get<std::string>("monitor", "monitor l2 update");
  if (monitor_name == "monitor residual") {
    monitor_type_ = Monitor::RESIDUAL;
    norm_type_ = MonitorNorm::ENORM;
  } else if (monitor_name == "monitor l2 residual") {
    monitor_type_ = Monitor::RESIDUAL;
    norm_type_ = MonitorNorm::L2;
  } else if (monitor_name == "monitor linf residual") {
    monitor_type_ = Monitor::RESIDUAL;
    norm_type_ = MonitorNorm::LINF;
  // } else if (monitor_name == "monitor preconditioned residual") {
  //   monitor_type_ = SOLVER_MONITOR_PCED_RESIDUAL;
  //   norm_type_ = MonitorNorm::ENORM;
  // } else if (monitor_name == "monitor preconditioned l2 residual") {
  //   monitor_type_ = SOLVER_MONITOR_PCED_RESIDUAL;
  //   norm_type_ = MonitorNorm::L2;
  // } else if (monitor_name == "monitor preconditioned linf residual") {
  //   monitor_type_ = SOLVER_MONITOR_PCED_RESIDUAL;
  //   norm_type_ = MonitorNorm::LINF;
  } else if (monitor_name == "monitor update") {
    monitor_type_ = Monitor::UPDATE;  // default value
    norm_type_ = MonitorNorm::ENORM;
  } else if (monitor_name == "monitor l2 update") {
    monitor_type_ = Monitor::UPDATE;  // default value
    norm_type_ = MonitorNorm::L2;
  } else if (monitor_name == "monitor linf update") {
    monitor_type_ = Monitor::UPDATE;  // default value
    norm_type_ = MonitorNorm::LINF;
  } else {
    Errors::Message m;
    m << "SolverNKA: Invalid monitor \"" << monitor_name << "\"";
    Exceptions::amanzi_throw(m);
  }
}


template <class Vector, class VectorSpace>
int
SolverDefault<Vector,VectorSpace>::Solve(const Teuchos::RCP<Vector>& u)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // reinitialize
  num_itrs_ = 0;
  pc_updates_ = 0;
  pc_calls_ = 0;
  function_calls_ = 0;
  residual_ = -1.0;
  status_ = MonitorStatus::CONTINUE;
  Restart_();

  // generate work space
  auto r = Teuchos::rcp(new Vector(u->getMap()));
  auto du = Teuchos::rcp(new Vector(u->getMap()));

  // variables to monitor the progress of the nonlinear solver
  error_ = std::numeric_limits<double>::max();
  residual_ = std::numeric_limits<double>::max();

  double previous_error(std::numeric_limits<double>::max());
  double l2_error(std::numeric_limits<double>::max());
  double l2_error_initial(std::numeric_limits<double>::max());
  double du_norm(std::numeric_limits<double>::max());
  double previous_du_norm(std::numeric_limits<double>::max());
  double du_norm_initial(std::numeric_limits<double>::max());
  double r_norm_initial(std::numeric_limits<double>::max());

  int divergence_count(0);
  int db_write_iter = 0;

  // start the iteration loop
  do {
    // Check for too many nonlinear iterations.
    if (num_itrs_ >= max_itrs_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Solve reached maximum of iterations (" << num_itrs_
                   << ")  error=" << error_ << " terminating..." << std::endl;
      status_ = MonitorStatus::MAX_ITERATIONS;
      return 1;
    }

    // Evaluate the nonlinear function.
    function_calls_++;
    fn_->Residual(u, r);
    residual_ = r->norm2();
    db_->WriteVector<Vector>(db_write_iter++, *r, u.ptr(), du.ptr());

    // If monitoring the residual, check for convergence.
    if (monitor_type_ == Monitor::RESIDUAL) {
      previous_error = error_;
      l2_error = residual_;

      if (norm_type_ == MonitorNorm::LINF)
        error_ = r->normInf();
      else if (norm_type_ == MonitorNorm::L2)
        error_ = residual_;
      else if (norm_type_ == MonitorNorm::ENORM)
        error_ = fn_->ErrorNorm(u, r);

      // We attempt to catch non-convergence early.
      // Stagnation based on L2 norm.
      if (num_itrs_ == 0) {
        l2_error_initial = l2_error;
      }

      // Check for convergence
      status_ = MonitorError_(error_, previous_error, l2_error, l2_error_initial, divergence_count);
      if (status_ == MonitorStatus::CONVERGED) return 0;
      else if (status_ != MonitorStatus::CONTINUE) return 1;
    }

    // Update the preconditioner if necessary.
    if (num_itrs_ % (pc_lag_ + 1) == 0) {
      pc_updates_++;
      fn_->UpdatePreconditioner(u);
    }

    // Apply the preconditioner to the nonlinear residual.
    pc_calls_++;
    int prec_error = fn_->ApplyPreconditioner(r, du);
    if (prec_error) {
      status_ = MonitorStatus::LINEAR_SOLVER_ERROR;
      return 1;
    }

    // Check for a method-specific direction.
    //
    // NOTE: It is up to the method to ensure that, after this call, if status
    // is CONTINUE, that du is valid, not diverging, and good to directly
    // apply.
    previous_du_norm = du_norm;
    std::tie(status_, du_norm) = ModifyCorrection_(r, u, du);
    if (status_ == MonitorStatus::CONVERGED) return 0;
    else if (status_ != MonitorStatus::CONTINUE) return 1;

    // Keep track of diverging iterations
    if (du_norm < 0) {
      // ModifyCorrection did not take the norm for us.
      du_norm = du->norm2();
    }

    // Next solution iterate and error estimate: u  = u - du
    u->update(-1.0, *du, 1.0);
    fn_->ChangedSolution();

    // Increment iteration counter.
    num_itrs_++;

    // If monitoring the update, check for convergence.
    if (monitor_type_ == Monitor::UPDATE) {
      previous_error = error_;
      l2_error = du_norm;

      if (norm_type_ == MonitorNorm::LINF)
        error_ = du->normInf();
      else if (norm_type_ == MonitorNorm::L2)
        error_ = du_norm;
      else if (norm_type_ == MonitorNorm::ENORM)
        error_ = fn_->ErrorNorm(u, du);

      // We attempt to catch non-convergence early.
      // Stagnation based on L2 norm.
      if (num_itrs_ == 1) {
        l2_error_initial = l2_error;
      }

      // Check for convergence
      status_ = MonitorError_(error_, previous_error, l2_error, l2_error_initial, divergence_count);
      if (status_ == MonitorStatus::CONVERGED) return 0;
      else if (status_ != MonitorStatus::CONTINUE) return 1;
    }
  } while(true); // continue do loop
}


template <class Vector, class VectorSpace>
MonitorStatus
SolverDefault<Vector,VectorSpace>::MonitorError_(double error,
        double previous_error, double l2_error, double l2_error_initial, int& divergence_count)
{
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << num_itrs_ << ": error=" << error << "  L2-error=" << l2_error << std::endl;

  if (error < tolerance_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Solver converged: " << num_itrs_
                 << " itrs, error=" << error << std::endl;
    return MonitorStatus::CONVERGED;

  } else if (num_itrs_ > max_itrs_stagnation_ && l2_error > l2_error_initial) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) 
      *vo_->os() << "Solver stagnating, L2-error=" << l2_error
                 << " > " << l2_error_initial << " (initial L2-error)" << std::endl;
    return MonitorStatus::STAGNATING;

  } else if (error > diverged_tol_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Solve failed, error " << error << " > " << diverged_tol_
                 << " (diverged)" << std::endl;
    return MonitorStatus::DIVERGED;

  } else if ((num_itrs_ > 1) &&
             (error > max_error_growth_factor_ * previous_error)) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Solver threatens to diverge, error " << error << " > "
                 << previous_error << " (previous error)" << std::endl;
    return MonitorStatus::DIVERGED;

  } else if (error >= previous_error) {
    divergence_count++;

    // If it does not recover quickly, abort.
    if (divergence_count >= max_divergence_count_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Solver is diverging repeatedly, terminating..." << std::endl;
      return MonitorStatus::DIVERGING;
    }

  } else {
    divergence_count = 0;
  }

  return MonitorStatus::CONTINUE;
}


} // namespace AmanziSolvers
} // namespace Amanzi

#endif
