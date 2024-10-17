/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

//! Adaptive timestep control based upon previous iteration count.
/*!

This is under development and is based on a posteriori error estimates.

*/

/*

DEVELOPER NOTE: This class has some issues -- it mixes between both
user-provided parameters for reduction factor and
DT_CONTROLLER_ADAPTIVE_REDUCTION, a hard-coded value.  It reads an
increase_factor_, but never uses it.  Unclear whether this is
working as intended or not... --ETC

*/

#ifndef AMANZI_ADAPTIVE_TIMESTEP_CONTROLLER_HH_
#define AMANZI_ADAPTIVE_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"

#include "TimestepController.hh"
#include "TimeIntegrationDefs.hh"

namespace Amanzi {

static const double DT_CONTROLLER_ADAPTIVE_INCREASE = 4.0;
static const double DT_CONTROLLER_ADAPTIVE_REDUCTION = 0.1;
static const double DT_CONTROLLER_ADAPTIVE_SAFETY_FACTOR = 0.9;
static const double DT_CONTROLLER_ADAPTIVE_ERROR_TOLERANCE = 1e-10;

template <class Vector>
class TimestepControllerAdaptive : public TimestepControllerRecoverable {
 public:
  TimestepControllerAdaptive(const std::string& name,
                             Teuchos::ParameterList& plist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<const Vector>& udot,
                             const Teuchos::RCP<const Vector>& udot_prev);

 protected:
  double getTimestep_(double dt, int iterations, bool valid) override;

 protected:
  int max_its_;
  int min_its_;
  double reduction_factor_, increase_factor_;

  Teuchos::RCP<const Vector> udot_prev_, udot_; // for error estimate
  double atol_, rtol_, p_;                // error parameters
};


inline double
getTimestepBase_(double dt,
                 const Epetra_MultiVector& u0,
                 const Epetra_MultiVector& u1,
                 double tol)
{
  double error_max = 0.0;
  double dTfactor(100.0), dTfactor_cell;

  int ncells_owned = u1.MyLength();
  for (int c = 0; c < ncells_owned; c++) {
    double error = fabs(u1[0][c] - u0[0][c]) * dt / 2;

    dTfactor_cell = sqrt(tol / std::max(error, DT_CONTROLLER_ADAPTIVE_ERROR_TOLERANCE));
    dTfactor = std::min(dTfactor, dTfactor_cell);

    error_max = std::max(error_max, error - tol);
  }

  dTfactor *= DT_CONTROLLER_ADAPTIVE_SAFETY_FACTOR;
  dTfactor = std::min(dTfactor, DT_CONTROLLER_ADAPTIVE_INCREASE);
  dTfactor = std::max(dTfactor, DT_CONTROLLER_ADAPTIVE_REDUCTION);

  double dT_tmp = dTfactor;
  auto comm = getCommWrapper(u0.Comm());
  comm->MinAll(&dT_tmp, &dTfactor, 1);

  double error_tmp = error_max;
  comm->MaxAll(&error_tmp, &error_max, 1);
  return dt * dTfactor;
}


/* ******************************************************************
* Constructor
****************************************************************** */
template <class Vector>
TimestepControllerAdaptive<Vector>::TimestepControllerAdaptive(const std::string& name,
        Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<const Vector>& udot,
        const Teuchos::RCP<const Vector>& udot_prev)
  : TimestepControllerRecoverable(name, plist, S),
    udot_prev_(udot_prev),
    udot_(udot)
{
  max_its_ = plist.get<int>("max iterations");
  min_its_ = plist.get<int>("min iterations");
  AMANZI_ASSERT(max_its_ > min_its_);
  AMANZI_ASSERT(min_its_ >= 0);

  reduction_factor_ = plist.get<double>("timestep reduction factor");
  AMANZI_ASSERT(reduction_factor_ >= 0.0);
  AMANZI_ASSERT(reduction_factor_ <= 1.0);

  increase_factor_ = plist.get<double>("timestep increase factor");
  AMANZI_ASSERT(increase_factor_ >= 1.0);

  // default value are fitted to atmospheric pressure.
  p_ = plist.get<double>("reference value", 101325.0);
  rtol_ = plist.get<double>("relative tolerance", 1.0e-4);
  atol_ = plist.get<double>("absolute tolerance", 10.0);
}


/* ******************************************************************
* Estimate new timestep by comparing the 1st and 2nd order time
* approximations.
****************************************************************** */
template <class Vector>
double
TimestepControllerAdaptive<Vector>::getTimestep_(double dt, int iterations, bool valid)
{
  if (iterations < 0 || iterations > max_its_ || !valid) {
    return dt * reduction_factor_;
  }
  double tol = rtol_ * p_ + atol_;
  return getTimestepBase_(dt, *udot_, *udot_prev_, tol);
}

template<>
inline double
TimestepControllerAdaptive<CompositeVector>::getTimestep_(double dt, int iterations, bool valid)
{
  if (iterations < 0 || iterations > max_its_ || !valid) {
    return dt * reduction_factor_;
  }
  double tol = rtol_ * p_ + atol_;
  return getTimestepBase_(dt, *udot_->ViewComponent("cell", false), *udot_prev_->ViewComponent("cell", false), tol);
}

template<>
inline double
TimestepControllerAdaptive<TreeVector>::getTimestep_(double dt, int iterations, bool valid)
{
  if (iterations < 0 || iterations > max_its_ || !valid) {
    return dt * reduction_factor_;
  }
  double tol = rtol_ * p_ + atol_;
  return getTimestepBase_(dt, *udot_->Data()->ViewComponent("cell", false), *udot_prev_->Data()->ViewComponent("cell", false), tol);
}

} // namespace Amanzi

#endif
