/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _NLScontrol_H_
#define _NLScontrol_H_

#include <MFTower.H>
#include <RSdata.H>
#include <NLSstatus.H>

struct RSdata; // Forward ref

struct NLScontrol
{
  /*
    Control structure for generic time-stepper algorithm with Newton-based nonlinear solver
   */
  NLScontrol();

  void SetMaxNewtonIterations(int max_iter);
  void SetMaxConsecutiveFails(int max_num);
  void SetDtRetryFactor(Real factor);
  void SetMaxConsecutiveFails2(int max_num);
  void SetDtRetryFactor2(Real factor);
  void SetDtRetryFactorF(Real factor);
  void SetMaxConsecutiveErrIncrease(int max_num);
  void SetConsecutiveErrIncreaseDtReduction(Real factor);

  void SetMaxNewtonIterationsForDt(int max_iter);
  void SetMaxConsecutiveSuccess(int max_num);
  void SetDtIncreaseFactor(Real factor);
  void SetMinNewtonIterationsForDt(int max_iter);
  void SetDtIncreaseFactor2(Real factor);
  void SetMinNewtonIterationsForDt2(int max_iter);
  void SetDtReductionFactor(Real factor);
  void SetMaxDt(Real dt_max);

  int  NLIterationsTaken() {return nl_iterations_taken;}
  void SetNLIterationsTaken(int iters);
  int  MaxNLIterations() {return max_nl_iterations;}
  void ResetCounters();

  bool AdjustDt(Real      dt,
		NLSstatus nl_solver_status,
		Real&     new_dt);




  /*
    Time-step adjustment strategy (dt_new given dt)
    -----------------------------------------------

    Initialize dt_new to dt
    if first:
       num_consecutive_increases = num_consecutive_success = 0, prev_abs_err = -1, first = false

    if nl solved:
       last_chance = FALSE

       if nl_iterations_taken < min_nl_iterations_for_dt:
         num_consecutive_success++, num_consecutive_increases = 0

         if num_consecutive_success >= max_num_consecutive_success:
           fac = time_step_increase_factor

           if nl_iterations_taken < min_nl_iterations_for_dt_2:
             fac = time_step_increase_factor_2

           dt_new = dt * fac

       if nl_iterations_taken > max_nl_iterations_for_dt:
         (force J recalc next try)
         num_consecutive_increases++, num_consecutive_success=0

         if nl_iterations_taken > max_nl_iterations_for_dt:
           dt_new = dt * time_step_reduction_factor

       num_consecutive_failures_1 = num_consecutive_failures_2 = 0

    else // step was rejected

      if num_consecutive_failures_1 <= max_num_consecutive_failures_1:
        dt_new = dt * time_step_retry_factor

      else:
        num_consecutive_failures_2++

        if num_consecutive_failures_2 <= max_num_consecutive_failures_2:
          dt_new = dt * time_step_retry_factor_2

        else:

          if last_chance:
            return FALSE

          dt_new = dt * time_step_retry_factor_f

          last_chance = TRUE

      num_consecutive_success = 0

      (force J recalc next try)

    dt_new = std::min(max_time_step_size,dt_new)

    return TRUE
   */

  // Richard solver data
  RSdata* rs_data;

  // Control parameters (see strategy above)
  std::string label;            // An identifier for this set of parameters
  int verbosity;                // Verbosity of this non-linear solver setup

  bool use_PETSc_snes;          // Use PETSc to solve the nonlinear problem
  bool use_fd_jac;              // Use matrix-free or finite-difference Jacobian
  bool use_dense_Jacobian;      // Use dense storage for Jacobian
  Real errfd;                   // Scale of perturbuation for fd Jacobian
  int  max_ls_iterations;       // Max line search attempts
  Real min_ls_factor;           // Minimum allowable line search scaling parameter
  Real ls_acceptance_factor;    // Maximun allowable reduction of residual in line search
  Real ls_reduction_factor;     // Scale to reduce line search scaling parameter for each try
  int  monitor_line_search;     // Write status information for line search
  int  maxit;                   // Maximum number of Newton iterations
  int  maxf;                    // Maximum number of function evaluations
  Real atol;                    // Absolute tolerance for norm(F)
  Real rtol;                    // Relative tolerance for norm(F)
  Real stol;                    // Absolute tolerance for norm(dp), dp is the Newton update
  bool scale_soln_before_solve; // Switch, scale P by 1/sigma before solve?
  bool centered_diff_J;         // Numerical J is computed with centered difference instead of forward
  bool ls_success;              // Pass-thru status from linesearch to convergence check
  std::string ls_reason;        // Status of linesearch to convergence success

  int max_nl_iterations;        // Max number of nonlinear (quasi-Newton iterations) to attempt
  Real max_nl_residual_norm;    // Largest workable residual value, abort solver if bigger
  int max_num_consecutive_success; // Number consecutive successful solves triggering increase dt
  int max_num_consecutive_failures_1; // Max allowable consecutive failures before retrying with dt reduced by time_step_retry_factor
  int max_num_consecutive_failures_2; // Max allowable consecutive failures before retrying with dt reduced by time_step_retry_factor_2
  int max_num_consecutive_increases; // Max allowable increases in residual beyond acceptance factor before retyring with dt reduced by consecutive_increase_reduction_factor
  Real consecutive_increase_reduction_factor;
  int min_nl_iterations_for_dt; // Maximum iterations to solve, below which next dt will be increased by time_step_increase_factor
  int min_nl_iterations_for_dt_2; // Maximum iterations to solve, below which next dt will be increased by time_step_increase_factor_2
  int max_nl_iterations_for_dt; // Maximum iterations to solve, above which next dt will not be modified
  Real max_time_step_size;      // Maximum time step to attempt
  Real time_step_increase_factor, time_step_increase_factor_2, time_step_reduction_factor;
  Real time_step_retry_factor, time_step_retry_factor_2, time_step_retry_factor_f;

private:

  // Status
  int nl_iterations_taken;
  Real nl_residual_norm;
  Real prev_abs_err;
  bool first;
  bool last_chance;
  int num_consecutive_success, num_consecutive_failures_1, num_consecutive_failures_2;
  int num_consecutive_increases;
};

std::ostream& operator<<(std::ostream& os, const NLScontrol& rhs);

#endif
