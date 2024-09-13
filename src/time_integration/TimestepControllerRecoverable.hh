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
#include "TimestepController.hh"

namespace Amanzi {

class State;

class TimestepControllerRecoverable : public TimestepController
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


} // namespace Amanzi


