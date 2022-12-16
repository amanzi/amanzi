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
#include "VerboseObject.hh"

namespace Amanzi {

TimeStepManager::TimeStepManager()
{
  dt_stable_storage = -1.;
  vo_ = Teuchos::rcp(new VerboseObject("TimeStepManager", Teuchos::ParameterList()));
}

TimeStepManager::TimeStepManager(Teuchos::ParameterList& verb_list)
{
  dt_stable_storage = -1.;
  vo_ = Teuchos::rcp(new VerboseObject("TimeStepManager", verb_list));
}

TimeStepManager::TimeStepManager(Teuchos::RCP<VerboseObject> vo_cd)
{
  dt_stable_storage = -1.;
  vo_ = vo_cd;
}

void
TimeStepManager::RegisterTimeEvent(double start, double period, double stop, bool phys)
{
  timeEvents_.push_back(TimeEvent(start, period, stop, phys));
}

void
TimeStepManager::RegisterTimeEvent(std::vector<double> times, bool phys)
{
  // make sure we only admit sorted arrays with unique entries
  std::vector<double> loc_times;
  loc_times = times;
  std::sort(loc_times.begin(), loc_times.end());
  loc_times.erase(std::unique(loc_times.begin(), loc_times.end(), near_equal), loc_times.end());
  timeEvents_.push_back(TimeEvent(loc_times, phys));
}

void
TimeStepManager::RegisterTimeEvent(double time, bool phys)
{
  timeEvents_.push_back(TimeEvent(time, phys));
}

double
TimeStepManager::TimeStep(double T, double dT, bool after_failure)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  Utils::Units units("molar");
  double next_T_all_events(1e99);
  bool physical = true;

  if (after_failure) dt_stable_storage = -1.;

  if ((dt_stable_storage > 0) && (!after_failure)) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Proposed dT=" << units.OutputTime(dT)
                 << ", stability from previous step changes it to "
                 << units.OutputTime(dt_stable_storage) << std::endl;
    }
    dT = dt_stable_storage;
    dt_stable_storage = -1.;
  }

  // loop over all events to find the next event time
  for (std::list<TimeEvent>::const_iterator i = timeEvents_.begin(); i != timeEvents_.end(); ++i) {
    double next_T_this_event(1e99);
    // there are two possible types of events
    switch (i->Type()) {
    case TimeEvent::SPS: {
      if (T < i->start_ && !Amanzi::near_equal(T, i->start_)) {
        next_T_this_event = i->start_;
      } else if ((i->stop_ == -1.0) || (T <= i->stop_ || Amanzi::near_equal(T, i->stop_))) {
        double n_periods = floor((T - i->start_) / i->period_);
        if (Amanzi::near_equal(n_periods + 1.0, (T - i->start_) / i->period_)) n_periods += 1.0;
        double tmp = i->start_ + (n_periods + 1.0) * i->period_;
        if (i->stop_ == -1.0 || tmp <= i->stop_) next_T_this_event = tmp;
      }
    }
    case TimeEvent::TIMES: {
      // we assume that i->times_ is ordered with unique elements
      for (std::vector<double>::const_iterator j = i->times_.begin(); j != i->times_.end(); ++j) {
        if (*j > T) {
          next_T_this_event = *j;
          break;
        }
      }
    }
    }
    //next_T_all_events = std::min(next_T_all_events, next_T_this_event);
    if (next_T_this_event < next_T_all_events) {
      physical = i->isPhysical();
      next_T_all_events = next_T_this_event;
    }
  }

  if (next_T_all_events == 1e99) return dT;
  double time_remaining(next_T_all_events - T);

  if (close(dT, time_remaining, 1.e-8)) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Proposed dT=" << units.OutputTime(dT)
                 << ", is near equal to next event time remaining "
                 << units.OutputTime(time_remaining) << "." << std::endl;
    }
    return time_remaining;

  } else if (dT > time_remaining) {
    if (!physical) dt_stable_storage = dT;
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Proposed dT=" << units.OutputTime(dT) << ", events limit it to "
                 << units.OutputTime(time_remaining) << std::endl;
    }
    return time_remaining;

  } else if (dT > 0.75 * time_remaining) {
    if (!physical) dt_stable_storage = dT;
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
  for (std::list<TimeEvent>::const_iterator i = timeEvents_.begin(); i != timeEvents_.end(); ++i) {
    switch (i->Type()) {
    case TimeEvent::SPS: {
      double time = i->start_;
      while (time <= end && (i->stop_ == -1.0 || time <= i->stop_)) {
        print_times.push_back(time);
        time += i->period_;
      }
    }
    case TimeEvent::TIMES: {
      for (std::vector<double>::const_iterator j = i->times_.begin(); j != i->times_.end(); ++j) {
        if (*j >= start && *j <= end) { print_times.push_back(*j); }
      }
    }
    }
  }
  std::sort(print_times.begin(), print_times.end());
  print_times.erase(std::unique(print_times.begin(), print_times.end(), near_equal),
                    print_times.end());
  for (std::vector<double>::const_iterator i = print_times.begin(); i != print_times.end(); ++i) {
    os << *i << " ";
  }
}

} // namespace Amanzi
