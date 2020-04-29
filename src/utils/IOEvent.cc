/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Markus Berndt
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "dbc.hh"
#include "errors.hh"
#include "TimeStepManager.hh"
#include "IOEvent.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Standard constructor.
// -----------------------------------------------------------------------------
IOEvent::IOEvent(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : plist_(plist), disabled_(false), units_()
{
  ReadParameters_();
};


// -----------------------------------------------------------------------------
// Constructor for a disabled Event.
// -----------------------------------------------------------------------------
IOEvent::IOEvent() : disabled_(true), units_(){};


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
IOEvent::RegisterWithTimeStepManager(TimeStepManager& tsm) const
{
  if (times_.size() != 0) { tsm.RegisterTimeEvent(times_.toVector(), false); }
  if (times_sps_.size() != 0) {
    for (Teuchos::Array<Teuchos::Array<double>>::const_iterator sps =
           times_sps_.begin();
         sps != times_sps_.end();
         ++sps) {
      tsm.RegisterTimeEvent((*sps)[0], (*sps)[1], (*sps)[2], false);
    }
  }
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
    if (cycles_.size() > 0) {
      for (Teuchos::Array<int>::const_iterator i = cycles_.begin();
           i != cycles_.end();
           ++i) {
        if (cycle == *i) { return true; }
      }
    }
    if (cycles_sps_.size() != 0) {
      for (Teuchos::Array<Teuchos::Array<int>>::const_iterator sps =
             cycles_sps_.begin();
           sps != cycles_sps_.end();
           ++sps) {
        if ((((*sps)[2] < 0) || (cycle <= (*sps)[2])) && ((*sps)[1] > 0)) {
          if ((*sps)[0] <= cycle) {
            int cycle_loc = cycle - (*sps)[0];
            if (cycle_loc % (*sps)[1] == 0) return true;
          }
        }
      }
    }
  }
  // if none of the conditions apply we do not write a visualization dump
  return false;
}

// -----------------------------------------------------------------------------
// Does vis need to dump the state by time?
// -----------------------------------------------------------------------------
bool
IOEvent::DumpRequested(double time) const
{
  if (!is_disabled()) {
    if (times_.size() > 0) {
      for (Teuchos::Array<double>::const_iterator i = times_.begin();
           i != times_.end();
           ++i) {
        if (Amanzi::near_equal(*i, time)) return true;
      }
    }
    if (times_sps_.size() != 0) {
      for (Teuchos::Array<Teuchos::Array<double>>::const_iterator sps =
             times_sps_.begin();
           sps != times_sps_.end();
           ++sps) {
        if (Amanzi::near_equal(time, (*sps)[0])) return true;
        if ((time > (*sps)[0]) &&
            (((*sps)[2] == -1.0) || (time <= (*sps)[2]))) {
          double n_periods = round((time - (*sps)[0]) / (*sps)[1]);
          double next_time = (*sps)[0] + n_periods * (*sps)[1];
          if (Amanzi::near_equal(time, next_time)) return true;
        }
      }
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
  if (plist_->isParameter("cycles start period stop")) {
    cycles_sps_.push_back(
      plist_->get<Teuchos::Array<int>>("cycles start period stop"));
  }

  bool done(false);
  int count(0);
  while (!done) {
    std::stringstream pname;
    pname << "cycles start period stop " << count;
    if (plist_->isParameter(pname.str())) {
      cycles_sps_.push_back(plist_->get<Teuchos::Array<int>>(pname.str()));
      count++;
    } else {
      done = true;
    }
  }

  if (plist_->isParameter("cycles")) {
    cycles_ = plist_->get<Teuchos::Array<int>>("cycles");
  }

  if (plist_->isParameter("times start period stop")) {
    auto times_sps =
      plist_->get<Teuchos::Array<double>>("times start period stop");

    // convert units
    auto my_units =
      plist_->get<std::string>("times start period stop units", "s");
    ValidUnitOrThrow_(my_units);
    bool success;
    for (auto& time : times_sps)
      time =
        time >= 0. ? units_.ConvertTime(time, my_units, "s", success) : -1.;
    AMANZI_ASSERT(success);

    times_sps_.push_back(times_sps);
  }

  done = false;
  count = 0;
  while (!done) {
    std::stringstream pname;
    pname << "times start period stop " << count;
    if (plist_->isParameter(pname.str())) {
      auto times_sps = plist_->get<Teuchos::Array<double>>(pname.str());

      // convert units
      pname << " units";
      auto my_units = plist_->get<std::string>(pname.str(), "s");
      ValidUnitOrThrow_(my_units);
      bool success;
      for (auto& time : times_sps)
        time =
          time >= 0. ? units_.ConvertTime(time, my_units, "s", success) : -1.;
      AMANZI_ASSERT(success);

      times_sps_.push_back(times_sps);
      count++;
    } else {
      done = true;
    }
  }

  if (plist_->isParameter("times")) {
    times_ = plist_->get<Teuchos::Array<double>>("times");

    // convert units
    auto my_units = plist_->get<std::string>("times units", "s");
    ValidUnitOrThrow_(my_units);
    bool success;
    for (auto& time : times_)
      time = units_.ConvertTime(time, my_units, "s", success);
    AMANZI_ASSERT(success);
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


} // namespace Amanzi
