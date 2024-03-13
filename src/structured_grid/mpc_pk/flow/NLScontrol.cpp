/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <NLScontrol.H>

// Solver parameter defaults
static bool use_PETSc_snes_DEF                        = true;
static bool use_fd_jac_DEF                            = true;
static bool use_dense_Jacobian_DEF                    = false;
static Real errfd_DEF                                 = 1.e-10;
static int max_ls_iterations_DEF                      = 10;
static Real min_ls_factor_DEF                         = 1.e-8;
static Real ls_acceptance_factor_DEF                  = 1.4;
static Real ls_reduction_factor_DEF                   = 0.1;
static int monitor_line_search_DEF                    = 0;
static int  maxit_DEF                                 = 30;
static int maxf_DEF                                   = 1e8;
static Real atol_DEF                                  = 1e-10;
static Real rtol_DEF                                  = 1e-20;
static Real stol_DEF                                  = 1e-12;
static bool scale_soln_before_solve_DEF               = true;
static bool centered_diff_J_DEF                       = false;
static bool ls_success_DEF                            = false;
static std::string ls_reason_DEF                      = "Invalid";
static int max_nl_iterations_DEF                      = 20;
static Real max_nl_residual_norm_DEF                  = -1;
static int max_num_consecutive_success_DEF            = 0;
static int max_num_consecutive_failures_1_DEF         = 3;
static int max_num_consecutive_failures_2_DEF         = 4;
static int max_num_consecutive_increases_DEF          = 15;
static Real consecutive_increase_reduction_factor_DEF = 0.4;
static int min_nl_iterations_for_dt_DEF               = 6;
static int min_nl_iterations_for_dt_2_DEF             = 3;
static int max_nl_iterations_for_dt_DEF               = 10;
static Real time_step_increase_factor_DEF             = 1.5;
static Real time_step_increase_factor_2_DEF           = 2.0;
static Real time_step_reduction_factor_DEF            = 0.8;
static Real time_step_retry_factor_DEF                = 0.5;
static Real time_step_retry_factor_2_DEF              = 0.1;
static Real time_step_retry_factor_f_DEF              = 0.01;
static Real max_time_step_size_DEF                    = 1.e10;
static int num_consecutive_success_DEF                = 0;
static int num_consecutive_failures_1_DEF             = 0;
static int num_consecutive_failures_2_DEF             = 0;
static int num_consecutive_increases_DEF              = 0;

Real norm0 = -1;
static Real max_residual_growth_factor = 1.e8; // FIXME: set this with rsparams
static Real min_dt = 1.e-2; // FIXME: set this with rsparams


void
NLScontrol::SetNLIterationsTaken(int iters) {
  nl_iterations_taken = iters;
}

NLScontrol::NLScontrol()
  : rs_data(0)
{
  use_PETSc_snes = use_PETSc_snes_DEF;

  // Set default values for all parameters
  use_fd_jac = use_fd_jac_DEF;
  centered_diff_J = centered_diff_J_DEF;
  use_dense_Jacobian = use_dense_Jacobian_DEF;
  errfd = errfd_DEF;
  max_ls_iterations = max_ls_iterations_DEF;
  min_ls_factor = min_ls_factor_DEF;
  ls_acceptance_factor = ls_acceptance_factor_DEF;
  ls_reduction_factor = ls_reduction_factor_DEF;
  monitor_line_search = monitor_line_search_DEF;
  maxit = maxit_DEF;
  maxf = maxf_DEF;
  atol = atol_DEF;
  rtol = rtol_DEF;
  stol = stol_DEF;
  scale_soln_before_solve = scale_soln_before_solve_DEF;
  ls_success = ls_success_DEF;
  ls_reason = ls_reason_DEF;
  max_nl_iterations = max_nl_iterations_DEF;
  max_nl_residual_norm = max_nl_residual_norm_DEF;
  max_num_consecutive_success = max_num_consecutive_success_DEF;
  max_num_consecutive_failures_1 = max_num_consecutive_failures_1_DEF;
  max_num_consecutive_failures_2 = max_num_consecutive_failures_2_DEF;
  max_num_consecutive_increases = max_num_consecutive_increases_DEF;
  consecutive_increase_reduction_factor = consecutive_increase_reduction_factor_DEF;
  min_nl_iterations_for_dt = min_nl_iterations_for_dt_DEF;
  min_nl_iterations_for_dt_2 = min_nl_iterations_for_dt_2_DEF;
  max_nl_iterations_for_dt = max_nl_iterations_for_dt_DEF;
  time_step_increase_factor = time_step_increase_factor_DEF;
  time_step_increase_factor_2 = time_step_increase_factor_2_DEF;
  time_step_reduction_factor = time_step_reduction_factor_DEF;
  time_step_retry_factor = time_step_retry_factor_DEF;
  time_step_retry_factor_2 = time_step_retry_factor_2_DEF;
  time_step_retry_factor_f = time_step_retry_factor_f_DEF;
  max_time_step_size = max_time_step_size_DEF;
  num_consecutive_failures_1 = num_consecutive_failures_1_DEF;
  num_consecutive_failures_2 = num_consecutive_failures_2_DEF;
  num_consecutive_increases = num_consecutive_increases_DEF;

  ResetCounters();
  first = true;
}

void NLScontrol::SetMaxConsecutiveFails(int max_num) {max_num_consecutive_failures_1=max_num;}
void NLScontrol::SetDtRetryFactor(Real factor) {time_step_retry_factor = factor;}
void NLScontrol::SetMaxConsecutiveFails2(int max_num) {max_num_consecutive_failures_2=max_num;}
void NLScontrol::SetDtRetryFactor2(Real factor) {time_step_retry_factor_2 = factor;}
void NLScontrol::SetDtRetryFactorF(Real factor) {time_step_retry_factor_f = factor;}
void NLScontrol::SetMaxConsecutiveErrIncrease(int max_incr) {max_num_consecutive_increases=max_incr;}
void NLScontrol::SetConsecutiveErrIncreaseDtReduction(Real redux) {consecutive_increase_reduction_factor=redux;}
void NLScontrol::SetMaxConsecutiveSuccess(int max_num) {max_num_consecutive_success=max_num;}
void NLScontrol::SetMaxNewtonIterations(int max_iter) {max_nl_iterations=max_iter;}
void NLScontrol::SetMaxNewtonIterationsForDt(int max_iter) {max_nl_iterations_for_dt=max_iter;}
void NLScontrol::SetMinNewtonIterationsForDt(int min_iter) {min_nl_iterations_for_dt=min_iter;}
void NLScontrol::SetMinNewtonIterationsForDt2(int min_iter) {min_nl_iterations_for_dt_2=min_iter;}
void NLScontrol::SetDtIncreaseFactor(Real factor) {time_step_increase_factor=factor;}
void NLScontrol::SetDtIncreaseFactor2(Real factor) {time_step_increase_factor_2=factor;}
void NLScontrol::SetDtReductionFactor(Real factor) {time_step_reduction_factor=factor;}
void NLScontrol::SetMaxDt(Real dt_max) {max_time_step_size=dt_max;}

void
NLScontrol::ResetCounters()
{
  nl_iterations_taken = 0;
  nl_residual_norm = -1;
  last_chance = false;;
  prev_abs_err = -1;
}

bool
NLScontrol::AdjustDt(Real      dt,
		     NLSstatus nl_solver_status,
		     Real&     dt_new) // Note: return bool for whether run should stop
{
  dt_new = dt;
  if (first) {
    num_consecutive_increases = 0;
    num_consecutive_success = max_num_consecutive_success; // Should increase immediately, unless after failure
    first = false;
    prev_abs_err = -1;
  }
  if (nl_solver_status == NLSstatus::NLS_SUCCESS)
  {
    last_chance = false;

    // "success" is when the error is reduced using small number of iters
    // In this case, increment counter for this event, reset "increase" counter
    // If this keeps happening, increase dt and reset the counter for these events
    //  (when we do, if the problem was  particularly easy, increase dt dramatically)
    if (nl_iterations_taken < min_nl_iterations_for_dt ) {

      num_consecutive_success++;
      num_consecutive_increases = 0;

      if (num_consecutive_success >= max_num_consecutive_success)
      {
        Real fac = time_step_increase_factor;
        if (nl_iterations_taken < min_nl_iterations_for_dt_2) {
          fac = time_step_increase_factor_2;
        }
        dt_new = dt * fac;
      }
    }

    // "increase" is when large number of iters
    // In this case, increment counter for this event, guarantee recalc of J,
    // and reset "success" counter
    // If this keeps happening, reduce dt and reset the counter for these events
    if (nl_iterations_taken > max_nl_iterations_for_dt  )
    {
      if (rs_data != 0) {
        rs_data->ResetJacobianCounter();
      }
      num_consecutive_increases++;
      num_consecutive_success = 0;

      if (nl_iterations_taken > max_nl_iterations_for_dt)
      {
        if (rs_data != 0) {
          rs_data->ResetJacobianCounter();
        }
        dt_new = dt * time_step_reduction_factor;
      }
    }

    num_consecutive_failures_1 = 0;
    num_consecutive_failures_2 = 0;

  }
  else {

    // step was rejected
    num_consecutive_failures_1++;

    if (num_consecutive_failures_1 <= max_num_consecutive_failures_1)
    {
      dt_new = dt * time_step_retry_factor;
    }
    else
    {
      num_consecutive_failures_2++;

      if (num_consecutive_failures_2 <= max_num_consecutive_failures_2)
      {
        dt_new = dt * time_step_retry_factor_2;
      }
      else
      {
        if (last_chance)  return false;
        dt_new = dt * time_step_retry_factor_f;
        last_chance = true;
      }
    }

    num_consecutive_success = 0;
    if (rs_data != 0) {
      rs_data->ResetJacobianCounter();
    }
  }

  dt_new = std::min(max_time_step_size,dt_new);
  return true;
}

std::ostream& operator<<(std::ostream& os, const NLScontrol& rhs)
{
  os << "NLScontrol configuration:\n";
  os << "{\n";
  os << "  label = " << rhs.label << "\n";
  os << "  verbosity = " << rhs.verbosity << "\n";
  os << "  use_PETSc_snes = " << rhs.use_PETSc_snes << "\n";
  os << "  use_fd_jac = " << rhs.use_fd_jac << "\n";
  os << "  use_dense_Jacobian = " << rhs.use_dense_Jacobian << "\n";
  os << "  errfd = " << rhs.errfd << "\n";
  os << "  max_ls_iterations = " << rhs.max_ls_iterations << "\n";
  os << "  min_ls_factor = " << rhs.min_ls_factor << "\n";
  os << "  ls_acceptance_factor = " << rhs.ls_acceptance_factor << "\n";
  os << "  ls_reduction_factor = " << rhs.ls_reduction_factor << "\n";
  os << "  monitor_line_search = " << rhs.monitor_line_search << "\n";
  os << "  maxit = " << rhs.maxit << "\n";
  os << "  maxf = " << rhs.maxf << "\n";
  os << "  atol = " << rhs.atol << "\n";
  os << "  rtol = " << rhs.rtol << "\n";
  os << "  stol = " << rhs.stol << "\n";
  os << "  scale_soln_before_solve = " << rhs.scale_soln_before_solve << "\n";
  os << "  centered_diff_J = " << rhs.centered_diff_J << "\n";
  os << "  max_nl_iterations = " << rhs.max_nl_iterations << "\n";
  os << "  max_nl_residual_norm = " << rhs.max_nl_residual_norm << "\n";
  os << "  max_num_consecutive_success = " << rhs.max_num_consecutive_success << "\n";
  os << "  max_num_consecutive_failures_1 = " << rhs.max_num_consecutive_failures_1 << "\n";
  os << "  max_num_consecutive_failures_2 = " << rhs.max_num_consecutive_failures_1 << "\n";
  os << "  max_num_consecutive_increases = " << rhs.max_num_consecutive_increases << "\n";
  os << "  consecutive_increase_reduction_factor = " << rhs.consecutive_increase_reduction_factor << "\n";
  os << "  min_nl_iterations_for_dt = " << rhs.min_nl_iterations_for_dt << "\n";
  os << "  min_nl_iterations_for_dt_2 = " << rhs.min_nl_iterations_for_dt_2 << "\n";
  os << "  max_nl_iterations_for_dt = " << rhs.max_nl_iterations_for_dt << "\n";
  os << "  time_step_increase_factor = " << rhs.time_step_increase_factor << "\n";
  os << "  time_step_increase_factor_2 = " << rhs.time_step_increase_factor_2 << "\n";
  os << "  time_step_reduction_factor = " << rhs.time_step_reduction_factor << "\n";
  os << "  time_step_retry_factor = " << rhs.time_step_retry_factor << "\n";
  os << "  time_step_retry_factor_2 = " << rhs.time_step_retry_factor_2 << "\n";
  os << "  time_step_retry_factor_f = " << rhs.time_step_retry_factor_f << "\n";
  os << "}";

  return os;
}
