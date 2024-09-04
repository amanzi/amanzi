/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Daniil Svyatskiy
*/

/*!

Manages time events, along with PK timestep selections, to determine the
step size taken by the simulation.

Note that most options provided to this class are only used in unit or
regression testing.

.. _time-step-manager-spec:
.. admonition:: time-step-manager-spec

   * `"constant timestep [s]`" ``[double]`` **optional** If provided, ignore the
     PK and all events and take a uniform, prescribed time step.
   * `"prescribed time history [s]`" ``Array(string)`` **optional** Array
     of timesteps to take, ignores the PK and all events.
   * `"prescribed time history file name`" ``[string]`` **optional** Path to an h5 file
     containing a list of dts, ignores the PK and all events.
   * `"prescribed time history header`" ``[string]`` **"timesteps"** Name of the
     dataset containing the dts.

*/

#pragma once

#include <ostream>
#include <vector>

#include "VerboseObject.hh"

namespace Amanzi {
namespace Utils {

template<typename T> class Event;

class TimeStepManager {
 public:
  TimeStepManager();
  TimeStepManager(Teuchos::ParameterList& plist);
  TimeStepManager(Teuchos::RCP<VerboseObject> vo_cd);
  void RegisterEvent(double start, double period, double stop, bool phys = true);
  void RegisterEvent(std::vector<double> times, bool phys = true);
  void RegisterEvent(double time, bool phys = true) {
    RegisterEvent(std::vector<double>{time});
  }
  void RegisterEvent(const Teuchos::RCP<const Event<double>>& te) {
    time_events_.emplace_back(te);
  }

  double TimeStep(const double T, const double dT, bool after_failure = false);
  void print(std::ostream& os, double start, double end) const;

 protected:
  std::vector<Teuchos::RCP<const Event<double>>> time_events_;
  bool manual_override_;
  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace Utils
} // namespace Amanzi

