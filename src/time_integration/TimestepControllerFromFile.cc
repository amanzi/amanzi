/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Time Integration

*/

// Timestep controller which loads a timestep history from file.

#include "dbc.hh"
#include "HDF5Reader.hh"
#include "TimestepControllerFromFile.hh"

namespace Amanzi {

TimestepControllerFromFile::TimestepControllerFromFile(Teuchos::ParameterList& plist)
  : TimestepController(plist), current_(0)
{
  std::string filename = plist.get<std::string>("file name");
  std::string header = plist.get<std::string>("timestep header", "timesteps");

  HDF5Reader reader(filename);
  reader.ReadData(header, dt_history_);
  dt_init_ = dt_history_[0];
}


// single method for timestep control
double
TimestepControllerFromFile::get_timestep(double dt, int iterations)
{
  double new_dt = -1.0;
  if (current_ == 0) {
    // the first request for a timestep is AFTER the first step is attempted
    // with the PK's initial timestep.  If this was successful, we should skip
    // the first step.  If it failed, we should give the first step.
    if (iterations < 0) {
      // failed, provide the first step
      new_dt = dt_history_[current_];
      current_++;
    } else {
      // successful, check to make sure the first simulation worked that way too
      AMANZI_ASSERT(std::abs(dt_history_[current_] - dt) < 1.e-6);
      current_++; // skip the 0th
      new_dt = dt_history_[current_];
      current_++;
    }
  } else {
    // iterations < 0 implies failed timestep
    if (iterations < 0) {
      Errors::TimeStepCrash m("TimestepController: prescribed time step size failed.");
      Exceptions::amanzi_throw(m);
    } else if (current_ < dt_history_.size()) {
      new_dt = dt_history_[current_];
      current_++;
    } else if (current_ == dt_history_.size()) {
      // the last step of the simulation still asks for a "next" dt prior to
      // finding out that the simulation is done.  Delay for a step and hope
      // the run ends.
      new_dt = dt_history_[current_ - 1];
      current_++;
    } else {
      Errors::Message m(
        "TimestepController: file contains insufficient number of timestep values.");
      Exceptions::amanzi_throw(m);
    }
  }
  return new_dt;
}

} // namespace Amanzi
