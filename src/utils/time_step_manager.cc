#include "time_step_manager.hh"
#include <iostream>
#include <cmath>

namespace Amanzi {

void TimeStepManager::RegisterTimeEvent(double start, double period, double stop) {
  timeEvents_.push_back(TimeEvent(start, period, stop));
}

void TimeStepManager::RegisterTimeEvent(const std::vector<double>& cycles) {
  timeEvents_.push_back(TimeEvent(cycles));
}

double TimeStepManager::TimeStep(const double T, const double dT) {
  double next_T_all_events(1e99);
  // loop over all events to find the next event time
  for (std::list<TimeEvent>::const_iterator i = timeEvents_.begin(); i != timeEvents_.end(); ++i) {
    double next_T_this_event(1e99);
    // there are two possible types of events
    switch (i->Type()) {
      case TimeEvent::SPS : {
        if ( T < i->start_ ) {
          next_T_this_event = i->start_;
        } else if ( (i->stop_ == -1.0) || ( T <= i->stop_ - i->period_ ) ) {
          double n_periods = floor( (T - i->start_ )/i->period_ );
          next_T_this_event = i->start_ + (n_periods+1.0)*i->period_;
        }
      }
      case TimeEvent::TIMES: {
        
      }
    }
    next_T_all_events = std::min(next_T_all_events, next_T_this_event);
  }
  double time_remaining(next_T_all_events - T);
  if (dT >= time_remaining) {
    return time_remaining;
  } else if ( dT > 0.75*time_remaining) {
    return 0.5*time_remaining;
  } else {
    return dT;
  } 
}

}
