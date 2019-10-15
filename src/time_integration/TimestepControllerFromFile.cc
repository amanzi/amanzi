/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

// Timestep controller which loads a timestep history from file.

#include "dbc.hh"
#include "HDF5Reader.hh"
#include "TimestepControllerFromFile.hh"

namespace Amanzi {

TimestepControllerFromFile::TimestepControllerFromFile(
  Teuchos::ParameterList& plist)
  : current_(0)
{
  std::string filename = plist.get<std::string>("file name");
  std::string header = plist.get<std::string>("timestep header", "timesteps");

  HDF5Reader reader(filename);
  reader.ReadData(header, dt_history_);
}


// single method for timestep control
double
TimestepControllerFromFile::get_timestep(double dt, int iterations)
{
  if (current_ == 0) {
    // the first request for a timestep is AFTER the first step is attempted
    // with the PK's initial timestep.  If this was successful, we should skip
    // the first step.  If it failed, we should give the first step.
    if (iterations < 0) {
      // failed, provide the first step
      current_++;
      return dt_history_[current_ - 1];
    } else {
      // successful, check to make sure the first simulation worked that way too
      AMANZI_ASSERT(std::abs(dt_history_[current_] - dt) < 1.e-6);
      current_ += 2;
      return dt_history_[current_ - 1];
    }
  } else {
    // iterations < 0 implies failed timestep
    if (iterations < 0) {
      Errors::Message m("TimestepController: time step crash");
      Exceptions::amanzi_throw(m);
    }
    if (current_ >= dt_history_.size()) { return -1.0; }

    current_++;
    return dt_history_[current_ - 1];
  }
}

} // namespace Amanzi
