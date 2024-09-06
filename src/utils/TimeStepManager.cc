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

#include "dbc.hh"
#include "TimeStepManager.hh"
#include "Units.hh"
#include "Event.hh"
#include "Reader.hh"
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
  if (plist.isParameter("prescribed timesteps [s]")) {
    manual_dts_ = plist.get<Teuchos::Array<double>>("prescribed timesteps [s]").toVector();
    manual_dts_i_ = 0;
    manual_override_ = true;

  } else if (plist.isParameter("prescribed timesteps file name")) {
    std::string filename = plist.get<std::string>("prescribed timesteps file name");
    std::string header = plist.get<std::string>("prescribed timesteps header", "timesteps");
    auto reader = createReader(filename);

    reader->read(header, manual_dts_);
    manual_dts_i_ = 0;
    manual_override_ = true;
  }

  vo_ = Teuchos::rcp(new VerboseObject("TimeStepManager", plist));
}


TimeStepManager::TimeStepManager(Teuchos::RCP<VerboseObject> vo_cd)
  : manual_override_(false),
    vo_(vo_cd)
{}


void
TimeStepManager::RegisterTimeEvent(double start, double period, double stop, bool phys)
{
  time_events_.emplace_back(Teuchos::rcp(new EventSPS<double>(start, period, stop)));
}


void
TimeStepManager::RegisterTimeEvent(const std::vector<double>& times, bool phys)
{
  time_events_.emplace_back(Teuchos::rcp(new EventList<double>(times)));
}

void
TimeStepManager::RegisterTimeEvent(double time, bool phys)
{
  RegisterTimeEvent(std::vector<double>{time});
}

void
TimeStepManager::RegisterTimeEvent(const Teuchos::RCP<const Event<double>>& te)
{
  time_events_.emplace_back(te);
}


double
TimeStepManager::TimeStep(double T, double dT, bool after_failure)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  Utils::Units units("molar");

  if (manual_override_) {
    // under manual override, the PK dt is prescribed
    if (after_failure) {
      Errors::Message msg("TimeStepManager: manually prescribed timestep failed.");
      Exceptions::amanzi_throw(msg);
    }
    if (manual_dts_i_ != manual_dts_.size()) {
      dT = manual_dts_[manual_dts_i_++];
    }
  }

  // loop over all events to find the next event time
  double next_T_all_events(std::numeric_limits<double>::max());
  for (auto i : time_events_) {
    double next_T_this_event = i->getNext(T);
    if (next_T_this_event >= 0. && next_T_this_event < next_T_all_events) {
      next_T_all_events = next_T_this_event;
    }
  }

  if (next_T_all_events == std::numeric_limits<double>::max()) return dT;

  double time_remaining(next_T_all_events - T);

  if (isNearEqual(dT, time_remaining, 1e4 * Event_EPS<double>::value)) {
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
    if (i->contains(start)) print_times.push_back(start);

    double time = i->getNext(start);
    while (time >= 0 && time < end) {
      print_times.push_back(time);
      time = i->getNext(time);
    }
  }

  std::sort(print_times.begin(), print_times.end());
  print_times.erase(std::unique(print_times.begin(), print_times.end(),
          [=](double a, double b) { return isNearEqual<double>(a,b); }),
                    print_times.end());
  for (auto t : print_times) os << t << " ";
}

} // namespace Utils
} // namespace Amanzi
