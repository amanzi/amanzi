/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "State.hh"
#include "TimestepControllerRecoverable.hh"

namespace Amanzi {

TimestepControllerRecoverable::TimestepControllerRecoverable(const std::string& name,
        Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S)
  : name_(Keys::cleanName(name, true)),
    dt_name_(Keys::cleanName(name, true) + "_dt_internal"),
    S_(S)
{
  double dt_init(-1.0);
  if (plist.isParameter("initial timestep [s]")) {
    dt_init = plist.get<double>("initial timestep [s]");
  } else {
    // by default, allow the driver's initial timestep
    dt_init = plist.get<double>("initial timestep", -1.0);
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


double
TimestepControllerRecoverable::getTimestep(double dt, int iterations, bool valid)
{
  // negative dt is used by some steady-state PKs in Amanzi, but a second step
  // is never taken.  Just don't error.
  if (dt < 0) return dt;

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
  return dt_recommended;
}


double
TimestepControllerRecoverable::getInitialTimestep()
{
  return *dt_internal_;
}


} // namespace Amanzi
