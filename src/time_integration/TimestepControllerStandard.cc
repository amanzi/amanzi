/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Time Integration

  Simple timestep control based upon previous iteration count.
*/

#include "errors.hh"
#include "dbc.hh"
#include "TimestepControllerStandard.hh"

namespace Amanzi {

TimestepControllerStandard::TimestepControllerStandard(const std::string& name,
                                                       Teuchos::ParameterList& plist,
                                                       const Teuchos::RCP<State>& S)
  : TimestepControllerRecoverable(name, plist, S)
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
}


double
TimestepControllerStandard::getTimestep_(double dt, int iterations, bool valid)
{
  double dt_next(dt);
  // iterations < 0 implies failed timestep
  if (iterations < 0 || iterations > max_its_ || !valid) {
    dt_next = dt * reduction_factor_;
  } else if (iterations < min_its_) {
    dt_next = dt * increase_factor_;
  }
  return dt_next;
}


} // namespace Amanzi
