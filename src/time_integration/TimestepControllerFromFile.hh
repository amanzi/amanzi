/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Timestep controller which loads a timestep history from file.
/*
  Time Integration

*/


/*!

This loads a timestep history from a file, then advances the step size with
those values.  This is mostly used for testing purposes, where we need to force
the same timestep history as previous runs to do regression testing.  Otherwise
even machine roundoff can eventually alter number of iterations enough to alter
the timestep history, resulting in solutions which are enough different to
cause doubt over their correctness.

.. _timestep-controller-from-file-spec:
.. admonition:: timestep-controller-from-file-spec

    * `"file name`" ``[string]`` Path to hdf5 file containing timestep information.
    * `"timestep header`" ``[string]`` Name of the dataset containing the history of timestep sizes.

*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "dbc.hh"
#include "Reader.hh"
#include "TimestepController.hh"

namespace Amanzi {

template<typename Vector>
class TimestepControllerFromFile : public TimestepController<Vector> {
 public:
  TimestepControllerFromFile(Teuchos::ParameterList& plist);

  // single method for timestep control
  double getTimestep(double dt, int iterations, bool valid) override;

 protected:
  std::vector<double> dt_history_;
  int current_;
};


template<typename Vector>
TimestepControllerFromFile<Vector>::TimestepControllerFromFile(Teuchos::ParameterList& plist)
  : TimestepController<Vector>(), current_(0)
{
  std::string filename = plist.get<std::string>("file name");
  std::string header = plist.get<std::string>("timestep header", "timesteps");

  auto reader = createReader(filename);
  reader->read(header, dt_history_);
  if (dt_history_.size() == 0) {
    Errors::Message m;
    m << "TimestepController: file \"" << filename << "\" timestep header has 0 times.";
    Exceptions::amanzi_throw(m);
  }
}


// single method for timestep control
template<typename Vector>
double
TimestepControllerFromFile<Vector>::getTimestep(double dt, int iterations, bool valid)
{
  double new_dt = -1.0;
  // iterations < 0 implies failed timestep
  if (iterations < 0 || !valid) {
    Errors::TimestepCrash m("TimestepController: prescribed timestep size failed.");
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
  return new_dt;
}

} // namespace Amanzi
