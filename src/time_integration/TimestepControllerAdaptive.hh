/* -------------------------------------------------------------------------
  AMANZI

  Author: Ethan Coon
  Author: Konstantin Lipnikov

  Adaptive timestep control based upon previous iteration count.
------------------------------------------------------------------------- */

#ifndef AMANZI_ADAPTIVE_TIMESTEP_CONTROLLER_HH_
#define AMANZI_ADAPTIVE_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TimestepController.hh"
#include "TimeIntegrationDefs.hh"

namespace Amanzi {

template<class Vector>
class TimestepControllerAdaptive : public TimestepController {

 public:
  TimestepControllerAdaptive(Teuchos::ParameterList& plist,
                             Teuchos::RCP<Vector> udot, Teuchos::RCP<Vector> udot_prev);

  // single method for timestep control
  double get_timestep_old(double dt, int iterations);
  double get_timestep(double dt, int iterations);

 protected:
  Teuchos::ParameterList plist_;

  int max_its_;
  int min_its_;
  double reduction_factor_;
  double increase_factor_;
  double max_dt_;
  double min_dt_;

 private:
  Teuchos::RCP<Vector> udot_prev_, udot_;  // for error estimate 
};


/* ******************************************************************
* Constructor 
****************************************************************** */
template<class Vector>
TimestepControllerAdaptive<Vector>::TimestepControllerAdaptive(
    Teuchos::ParameterList& plist,
    Teuchos::RCP<Vector> udot, Teuchos::RCP<Vector> udot_prev)
    : plist_(plist), udot_(udot), udot_prev_(udot_prev)
{
  max_its_ = plist_.get<int>("max iterations");
  min_its_ = plist_.get<int>("min iterations");
  ASSERT(max_its_ > min_its_);
  ASSERT(min_its_ >= 0);

  reduction_factor_ = plist_.get<double>("time step reduction factor");
  ASSERT(reduction_factor_ >= 0.0);
  ASSERT(reduction_factor_ <= 1.0);

  increase_factor_ = plist_.get<double>("time step increase factor");
  ASSERT(increase_factor_ >= 1.0);

  max_dt_ = plist_.get<double>("max time step");
  min_dt_ = plist_.get<double>("min time step");
}


/* ******************************************************************
* Estimate new time step.
****************************************************************** */
template<class Vector>
double TimestepControllerAdaptive<Vector>::get_timestep_old(double dt, int iterations)
{
  double dt_next(dt);

  // iterations < 0 implies failed timestep
  if (iterations < 0 || iterations > max_its_) {
    dt_next = dt * reduction_factor_;
  } else if (iterations < min_its_) {
    dt_next = dt * increase_factor_;
  }

  // check max step size
  if (dt_next > max_dt_) dt_next = max_dt_;

  // check min step size
  if (dt_next < min_dt_) {
    if (iterations < 0) {
      std::string msg = "dT is too small, terminating...";
      Errors::Message m(msg);
      Exceptions::amanzi_throw(m);
    } else {
      dt_next = min_dt_;
    }
  }

  // check that if we have failed, our step size has decreased.
  if (iterations < 0) {
    if (dt - dt_next < 1.0e-10) {
      std::string msg = "dT change is too small, terminating...";
      Errors::Message m(msg);
      Exceptions::amanzi_throw(m);
    }
  }

  return dt_next;
}


/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations. 
****************************************************************** */
template<class Vector>
double TimestepControllerAdaptive<Vector>::get_timestep(double dt, int iterations)
{
  if (iterations < 0 || iterations > max_its_) {
    return dt * reduction_factor_;
  }
    
  double rtol(1.0), atol(1e-5), p(101325.0);
  double tol, error, error_max = 0.0;
  double dTfactor(100.0), dTfactor_cell;

#ifdef VECTORisEPETRA_VECTOR
  Epetra_MultiVector& u1 = *udot_;
  Epetra_MultiVector& u0 = *udot_prev_;
#else
  Epetra_MultiVector& u1 = *udot_->ViewComponent("cell");
  Epetra_MultiVector& u0 = *udot_prev_->ViewComponent("cell");
#endif

  int ncells_owned = u1.MyLength();

  for (int c = 0; c < ncells_owned; c++) {
    error = fabs(u1[0][c] - u0[0][c]) * dt / 2;
    tol = rtol * p + atol;

    dTfactor_cell = sqrt(tol / std::max(error, DT_CONTROLLER_ADAPTIVE_ERROR_TOLERANCE));
    dTfactor = std::min(dTfactor, dTfactor_cell);

    error_max = std::max(error_max, error - tol);
  }

  dTfactor *= DT_CONTROLLER_ADAPTIVE_SAFETY_FACTOR;
  dTfactor = std::min(dTfactor, DT_CONTROLLER_ADAPTIVE_INCREASE);
  dTfactor = std::max(dTfactor, DT_CONTROLLER_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
  double dT_tmp = dTfactor;
  udot_->Comm().MinAll(&dT_tmp, &dTfactor, 1);  // find the global minimum
 
  double error_tmp = error_max;
  udot_->Comm().MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif
  return dt * dTfactor;
}


/*
  // Calculate time derivative and 2nd-order solution approximation.
  // Estimate of a time step multiplier overrides the above estimate.
  if (ti_specs->dT_method == DT_CONTROLLER_ADAPTIVE) {
    const Epetra_MultiVector& pressure = *S_->GetFieldData("pressure")->ViewComponent("face");
    Epetra_MultiVector& p = *solution->ViewComponent("face");

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = (p[0][c] - pressure[0][c]) / dT; 
      p[0][c] = pressure[0][c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dT / 2;
    }

    double err, dTfactor;
    err = AdaptiveTimeStepEstimate(&dTfactor);
    if (err > 0.0) throw 1000;  // fix (lipnikov@lan.gov)
    dTnext = std::min(dT_MPC * dTfactor, ti_specs->dTmax);
  }
*/

} // namespace Amanzi

#endif
