/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Utils

  IO event -- base class for reading or writing data.  Mostly just manages when
  to do the I/O.
*/

#include "dbc.hh"
#include "errors.hh"
#include "TimeStepManager.hh"
#include "Event.hh"
#include "IOEvent.hh"

namespace Amanzi {
namespace Utils {

// -----------------------------------------------------------------------------
// Standard constructor.
// -----------------------------------------------------------------------------
IOEvent::IOEvent(Teuchos::ParameterList& plist) : plist_(plist), units_(), disabled_(false)
{
  ReadParameters_();
};


// -----------------------------------------------------------------------------
// Constructor for a disabled Event.
// -----------------------------------------------------------------------------
IOEvent::IOEvent() : units_(), disabled_(true){};


// -----------------------------------------------------------------------------
// Set flag to not act.
// -----------------------------------------------------------------------------
bool
IOEvent::is_disabled() const
{
  return disabled_;
}
void
IOEvent::disable(bool disabled)
{
  disabled_ = disabled;
}


// -----------------------------------------------------------------------------
// Place vis times in the manager.
// -----------------------------------------------------------------------------
void
IOEvent::RegisterWithTimeStepManager(TimeStepManager& tsm)
{
  for (const auto& te : time_events_) tsm.RegisterTimeEvent(te);
}


// -----------------------------------------------------------------------------
// Does vis need to dump the state?
// -----------------------------------------------------------------------------
bool
IOEvent::DumpRequested(int cycle, double time) const
{
  return DumpRequested(cycle) || DumpRequested(time);
}


// -----------------------------------------------------------------------------
// Does vis need to dump the state by cycle?
// -----------------------------------------------------------------------------
bool
IOEvent::DumpRequested(int cycle) const
{
  if (!is_disabled()) {
    for (const auto& ce : cycle_events_) {
      if (ce->contains(cycle)) return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
// Does vis need to dump the state by time?
// -----------------------------------------------------------------------------
bool
IOEvent::DumpRequested(double time) const
{
  if (!is_disabled()) {
    for (const auto& te : time_events_) {
      if (te->contains(time)) return true;
    }
  }
  return false;
}


// -----------------------------------------------------------------------------
// Parse the parameter list and init.
// -----------------------------------------------------------------------------
void
IOEvent::ReadParameters_()
{
  // cycle-based events
  if (plist_.isParameter("cycles start period stop")) {
    auto c_sps = plist_.get<Teuchos::Array<int>>("cycles start period stop");
    if (c_sps.size() != 3) {
      Errors::Message msg("Array of incorrect length provided to IOEvent \"cycles start period stop\"");
      Exceptions::amanzi_throw(msg);
    }
    cycle_events_.push_back(Teuchos::rcp(new EventSPS<int>(c_sps[0], c_sps[1], c_sps[2])));
  }

  bool done(false);
  int count(0);
  while (!done) {
    std::stringstream pname;
    pname << "cycles start period stop " << count;
    if (plist_.isParameter(pname.str())) {
      auto c_sps = plist_.get<Teuchos::Array<int>>(pname.str());
      if (c_sps.size() != 3) {
        Errors::Message msg;
        msg << "Array of incorrect length provided to IOEvent \""
            << pname.str() << "\"";
        Exceptions::amanzi_throw(msg);
      }
      cycle_events_.push_back(Teuchos::rcp(new EventSPS<int>(c_sps[0], c_sps[1], c_sps[2])));
      count++;
    } else {
      done = true;
    }
  }

  if (plist_.isParameter("cycles")) {
    auto cycles = plist_.get<Teuchos::Array<int>>("cycles");
    cycle_events_.push_back(Teuchos::rcp(new EventList<int>(cycles.toVector())));
  }

  // time-based events
  if (plist_.isParameter("times start period stop")) {
    auto t_sps = plist_.get<Teuchos::Array<double>>("times start period stop");
    if (t_sps.size() != 3) {
      Errors::Message msg("Array of incorrect length provided to IOEvent \"times start period stop\"");
      Exceptions::amanzi_throw(msg);
    }

    // convert units
    auto my_units = plist_.get<std::string>("times start period stop units", "s");
    ValidUnitOrThrow_(my_units);
    bool success;
    for (auto& time : t_sps)
      time = time >= 0. ? units_.ConvertTime(time, my_units, "s", success) : -1.;
    AMANZI_ASSERT(success);

    time_events_.push_back(Teuchos::rcp(new EventSPS<double>(t_sps[0], t_sps[1], t_sps[2])));
  }

  done = false;
  count = 0;
  while (!done) {
    std::stringstream pname;
    pname << "times start period stop " << count;
    if (plist_.isParameter(pname.str())) {
      auto t_sps = plist_.get<Teuchos::Array<double>>(pname.str());

      // convert units
      pname << " units";
      auto my_units = plist_.get<std::string>(pname.str(), "s");
      ValidUnitOrThrow_(my_units);
      bool success;
      for (auto& time : t_sps)
        time = time >= 0. ? units_.ConvertTime(time, my_units, "s", success) : -1.;
      AMANZI_ASSERT(success);

      time_events_.push_back(Teuchos::rcp(new EventSPS<double>(t_sps[0], t_sps[1], t_sps[2])));
      count++;
    } else {
      done = true;
    }
  }

  if (plist_.isParameter("times")) {
    auto times = plist_.get<Teuchos::Array<double>>("times");

    // convert units
    auto my_units = plist_.get<std::string>("times units", "s");
    ValidUnitOrThrow_(my_units);
    bool success;
    for (auto& time : times) time = units_.ConvertTime(time, my_units, "s", success);
    AMANZI_ASSERT(success);

    time_events_.push_back(Teuchos::rcp(new EventList<double>(times.toVector())));
  }
}


// -----------------------------------------------------------------------------
// Check a user-provided unit string and throw if it is not a valid time unit.
// -----------------------------------------------------------------------------
void
IOEvent::ValidUnitOrThrow_(const std::string& my_units)
{
  if (!units_.IsValidTime(my_units)) {
    Errors::Message msg;
    msg << "IOEvent: unknown time units type: \"" << my_units
        << "\"  Valid are: " << units_.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
}


} // namespace Utils
} // namespace Amanzi
