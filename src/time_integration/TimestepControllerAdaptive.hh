/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/
//! Adaptive timestep control based upon previous iteration count.

/*!

This is under development and is based on a posteriori error estimates.

*/

#ifndef AMANZI_ADAPTIVE_TIMESTEP_CONTROLLER_HH_
#define AMANZI_ADAPTIVE_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"

#include "TimestepController.hh"
#include "TimeIntegrationDefs.hh"

namespace Amanzi {

template <class Vector>
class TimestepControllerAdaptive : public TimestepController {
 public:
  TimestepControllerAdaptive(Teuchos::ParameterList& plist,
                             Teuchos::RCP<Vector> udot,
                             Teuchos::RCP<Vector> udot_prev);

  // single method for timestep control
  double get_timestep(double dt, int iterations);

 private:
  double get_timestep_base_(double dt, const Epetra_MultiVector& u0, const Epetra_MultiVector& u1);

 protected:
  Teuchos::ParameterList plist_;

  int max_its_;
  int min_its_;
  double reduction_factor_;
  double increase_factor_;
  double max_dt_;
  double min_dt_;

 private:
  Teuchos::RCP<Vector> udot_prev_, udot_; // for error estimate
  double atol_, rtol_, p_;                // error parameters
};


/* ******************************************************************
* Constructor 
****************************************************************** */
template <class Vector>
TimestepControllerAdaptive<Vector>::TimestepControllerAdaptive(Teuchos::ParameterList& plist,
                                                               Teuchos::RCP<Vector> udot,
                                                               Teuchos::RCP<Vector> udot_prev)
  : TimestepController(plist), plist_(plist), udot_prev_(udot_prev), udot_(udot)
{
  max_its_ = plist_.get<int>("max iterations");
  min_its_ = plist_.get<int>("min iterations");
  AMANZI_ASSERT(max_its_ > min_its_);
  AMANZI_ASSERT(min_its_ >= 0);

  reduction_factor_ = plist_.get<double>("time step reduction factor");
  AMANZI_ASSERT(reduction_factor_ >= 0.0);
  AMANZI_ASSERT(reduction_factor_ <= 1.0);

  increase_factor_ = plist_.get<double>("time step increase factor");
  AMANZI_ASSERT(increase_factor_ >= 1.0);

  if (plist_.isParameter("max time step [s]")) {
    max_dt_ = plist_.get<double>("max time step [s]");
  } else {
    max_dt_ = plist_.get<double>("max time step");
  }
  if (plist_.isParameter("min time step [s]")) {
    min_dt_ = plist_.get<double>("min time step [s]");
  } else {
    min_dt_ = plist_.get<double>("min time step");
  }

  // default value are fitted to atmospheric pressure.
  p_ = plist_.get<double>("reference value", 101325.0);
  rtol_ = plist_.get<double>("relative tolerance", 1.0e-4);
  atol_ = plist_.get<double>("absolute tolerance", 10.0);
}


/* ******************************************************************
* Estimate new time step by comparing the 1st and 2nd order time 
* approximations. 
****************************************************************** */
template <class Vector>
double
TimestepControllerAdaptive<Vector>::get_timestep(double dt, int iterations)
{
  std::string msg("TimestepControllerAdaptive is implemented for a limited set of vectors.");
  Errors::Message m(msg);
  Exceptions::amanzi_throw(m);
  return 0.0;
}


template <>
inline double
TimestepControllerAdaptive<Epetra_MultiVector>::get_timestep(double dt, int iterations)
{
  if (iterations < 0 || iterations > max_its_) { return dt * reduction_factor_; }

  Epetra_MultiVector& u1 = *udot_;
  Epetra_MultiVector& u0 = *udot_prev_;
  return get_timestep_base_(dt, u0, u1);
}


template <>
inline double
TimestepControllerAdaptive<CompositeVector>::get_timestep(double dt, int iterations)
{
  if (iterations < 0 || iterations > max_its_) { return dt * reduction_factor_; }

  Epetra_MultiVector& u1 = *udot_->ViewComponent("cell");
  Epetra_MultiVector& u0 = *udot_prev_->ViewComponent("cell");
  return get_timestep_base_(dt, u0, u1);
}


template <class Vector>
double
TimestepControllerAdaptive<Vector>::get_timestep_base_(double dt,
                                                       const Epetra_MultiVector& u0,
                                                       const Epetra_MultiVector& u1)
{
  double tol, error, error_max = 0.0;
  double dTfactor(100.0), dTfactor_cell;

  int ncells_owned = u1.MyLength();
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs(u1[0][c] - u0[0][c]) * dt / 2;
    tol = rtol_ * p_ + atol_;

    dTfactor_cell = sqrt(tol / std::max(error, DT_CONTROLLER_ADAPTIVE_ERROR_TOLERANCE));
    dTfactor = std::min(dTfactor, dTfactor_cell);

    error_max = std::max(error_max, error - tol);
  }

  dTfactor *= DT_CONTROLLER_ADAPTIVE_SAFETY_FACTOR;
  dTfactor = std::min(dTfactor, DT_CONTROLLER_ADAPTIVE_INCREASE);
  dTfactor = std::max(dTfactor, DT_CONTROLLER_ADAPTIVE_REDUCTION);

  double dT_tmp = dTfactor;
  auto comm = getCommWrapper(udot_->Comm());
  comm->MinAll(&dT_tmp, &dTfactor, 1);

  double error_tmp = error_max;
  comm->MaxAll(&error_tmp, &error_max, 1);

  return std::min(dt * dTfactor, max_dt_);
}

} // namespace Amanzi

#endif
