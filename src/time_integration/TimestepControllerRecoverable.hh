/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

A helper class for a TimestepController that can recover from event-based
timestep reductions.

This saves the internal dt in State, and can recover the preferred dt from
arbitrarily small steps introduced by events.

It also provides some basic functionality for adaptive timesteppers, setting a
max, min, and initial dt.

.. _timestep-controller-recoverable-spec
.. admonition:: timestep-controller-recoverable-spec

   * `"initial timestep [s]`" ``[double]`` The max timestep size allowed.
   * `"max timestep [s]`" ``[double]`` The max timestep size allowed.
   * `"min timestep [s]`" ``[double]`` The min timestep size allowed.  If the step
     has failed and the new step is below this cutoff, the simulation fails.

*/

#pragma once

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "State.hh"
#include "TimestepController.hh"

namespace Amanzi {

template <typename Vector>
class TimestepControllerRecoverable : public TimestepController<Vector>
{
 public:
  TimestepControllerRecoverable(const std::string& name,
          Teuchos::ParameterList& plist,
          const Teuchos::RCP<State>& S);

  double getTimestep(double dt, int iterations, bool valid) override final;
  double getInitialTimestep() override;

 protected:
  // client classes should implement this instead
  virtual double getTimestep_(double dt, int iterations, bool valid) = 0;

 protected:
  std::string name_, dt_name_;
  Teuchos::RCP<State> S_;
  Teuchos::RCP<double> dt_internal_;
  double dt_min_;
  double dt_max_;
};


template <typename Vector>
TimestepControllerRecoverable<Vector>::TimestepControllerRecoverable(const std::string& name,
        Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S)
  : TimestepController<Vector>(),
    name_(Keys::cleanName(name, true)),
    S_(S),
    dt_name_(Keys::cleanName(name, true) + "_dt_internal")
{
  double dt_init(-1.0);
  if (plist.isParameter("initial timestep [s]")) {
    dt_init = plist.get<double>("initial timestep [s]");
  } else {
    dt_init = plist.get<double>("initial timestep", 1.0);
  }

  if (plist.isParameter("max timestep")) {
    dt_max_ = plist.get<double>("max timestep");
  } else {
    dt_max_ = plist.get<double>("max timestep [s]", std::numeric_limits<double>::max());
  }

  if (plist.isParameter("min timestep [s]")) {
    dt_min_ = plist.get<double>("min timestep [s]");
  } else {
    dt_min_ = plist.get<double>("min timestep", 10*std::numeric_limits<double>::min());
  }

  // create state memory for internal dt and give it to state
  dt_internal_ = Teuchos::rcp(new double(dt_init));
  if (S_.get()) {
    S_->Require<double>(dt_name_, Tags::DEFAULT, name_);
    S_->SetPtr<double>(dt_name_, Tags::DEFAULT, name_, dt_internal_);
    S_->GetRecordW(dt_name_, Tags::DEFAULT, name_).set_initialized();
    S_->GetRecordW(dt_name_, Tags::DEFAULT, name_).set_io_checkpoint();
    S_->GetRecordW(dt_name_, Tags::DEFAULT, name_).set_io_vis(false);
  }
}


template <typename Vector>
double
TimestepControllerRecoverable<Vector>::getTimestep(double dt, int iterations, bool valid)
{
  double dt_prev_recommended = *dt_internal_;
  double dt_recommended = getTimestep_(dt, iterations, valid);

  if (dt <= dt_prev_recommended &&
      dt <= dt_recommended &&
      dt_recommended < dt_prev_recommended) {
    // We took a smaller step than we recommended, likely due to
    // constraints from other PKs or events like vis (dt <= dt_internal),
    // and it worked well enough that the newly recommended step size
    // didn't decrease (dt <= dt_solver).  Do not reduce our
    // recommendation.
    dt_recommended = dt_prev_recommended;
  }

  dt_recommended = std::min(dt_max_, dt_recommended);

  if (dt_recommended < dt_min_) {
    Errors::TimestepCrash msg;
    msg << "TimestepController reduced the timestep to " << dt_recommended << ", which is less than the minimum allowed dt, " << dt_min_;
    Exceptions::amanzi_throw(msg);
  }

  *dt_internal_ = dt_recommended;
  std::cout << "getTimestep: dt = " << dt << ", dt_prev_recommended = " << dt_prev_recommended << ", dt_recommended = " << dt_recommended << std::endl;
  return dt_recommended;
}


template <typename Vector>
double
TimestepControllerRecoverable<Vector>::getInitialTimestep()
{
  return *dt_internal_;
}


} // namespace Amanzi


