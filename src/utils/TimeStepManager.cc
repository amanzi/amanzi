/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Daniil Svyatskiy
*/

/*
  Utils

*/

#include <algorithm>
#include <iostream>

#include "TimeStepManager.hh"
#include "Units.hh"
#include "Event.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Utils {

TimeStepManager::TimeStepManager()
  : manual_override_(false)
{
  vo_ = Teuchos::rcp(new VerboseObject("TimeStepManager", Teuchos::ParameterList()));
}

TimeStepManager::TimeStepManager(Teuchos::ParameterList& plist)
  : manual_override_(false)
{
  // manual override
  if (plist.hasParameter("constant timestep [s]")) {
    double dt = plist.get<double>("constant timestep [s]");
    double t0 = plist.get<double>("constant timestep start time [s]", 0.0);
    RegisterTimeEvent(t0, dt, -1., false);
    manual_override_ = true;

  } else if (plist.hasParameter("prescribed time history [s]")) {
    auto dt_manual = plist.get<Teuchos::Array<std::string>>("prescribed time history [s]");
    RegisterTimeEvent(dt_manual.toVector(), false);
    manual_override_ = true;

  } else if (plist.hasParameter("prescribed time file name")) {
    std::string fname = plist.get<std::string>("prescribed time history file name");
    std::string dset = plist.get<std::string>("prescribed time history header", "timesteps");
    auto reader = createReader(filename);

    Teuchos::Array<double> dts;
    reader->read(header, dts);
    RegisterTimeEvent(dts.toVector(), false);
    manual_override_ = true;
  }

  vo_ = Teuchos::rcp(new VerboseObject("TimeStepManager", verb_list));
}


TimeStepManager::TimeStepManager(Teuchos::RCP<VerboseObject> vo_cd)
  : manual_override_(false),
    vo_(vo_cd)
{}


void
TimeStepManager::RegisterTimeEvent(double start, double period, double stop, bool phys)
{
  if (!manual_override_)
    timeEvents_.emplace_back(Teuchos::rcp(new TimeEvent<double>(start, period, stop)));
}


void
TimeStepManager::RegisterTimeEvent(std::vector<double> times, bool phys)
{
  if (!manual_override_)
    timeEvents_.emplace_back(Teuchos::rcp(new TimeEvent<double>(times)));
}


void
TimeStepManager::RegisterTimeEvent(double time, bool phys)
{
  if (!manual_override_)
    timeEvents_.emplace_back(Teuchos::rcp(new TimeEvent<double>(time, phys)));
}


double
TimeStepManager::TimeStep(double T, double dT, bool after_failure)
{
  if (manual_override_) {
    if (after_failure) {
      Errors::Message msg("Manually override timestep failed.");
      Exceptions::amanzi_throw(msg);
    }
    AMANZI_ASSERT(time_events_.size() == 1);
    return time_events_[0].getNext(T);
  }

  Teuchos::OSTab tab = vo_->getOSTab();
  Utils::Units units("molar");
  double next_T_all_events(std::numeric_limits<double>::max;);

  // loop over all events to find the next event time
  for (auto i : time_events_) {
    double next_T_this_event = i->getNext(T);
    if (next_T_this_event >= 0. && next_T_this_event < next_T_all_events) {
      next_T_all_events = next_T_this_event;
    }
  }

  if (next_T_all_events == std::numeric_limits<double>::max) return dT;
  double time_remaining(next_T_all_events - T);

  if (isNearEqual(dT, time_remaining, 1e4*TSM_EPS)) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Proposed dT=" << units.OutputTime(dT)
                 << ", is near equal to next event time remaining "
                 << units.OutputTime(time_remaining) << "." << std::endl;
    }
    return time_remaining;

  } else if (dT > time_remaining) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Proposed dT=" << units.OutputTime(dT) << ", events limit it to "
                 << units.OutputTime(time_remaining) << std::endl;
    }
    return time_remaining;

  } else if (dT > 0.75 * time_remaining) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Proposed dT=" << units.OutputTime(dT) << ", events limit it to "
                 << units.OutputTime(0.5 * time_remaining) << std::endl;
    }
    return 0.5 * time_remaining;

  } else {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Accepted proposed dT=" << units.OutputTime(dT) << std::endl;
    }
    return dT;
  }
}


void
TimeStepManager::print(std::ostream& os, double start, double end) const
{
  // create a sorted array of the times between start and end and print it
  std::vector<double> print_times;
  for (auto i : time_events_) {
    double time = start;
    while (time > 0 && time < end) {
      time = i->getNext(time);
      print_times.push_back(time);
    }
  }

  std::sort(print_times.begin(), print_times.end());
  print_times.erase(std::unique(print_times.begin(), print_times.end(), isNearEqual),
                    print_times.end());
  for (auto t : print_times) os << t << " ";
}

} // namespace Utils
} // namespace Amanzi
