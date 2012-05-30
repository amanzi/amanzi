#include "time_step_manager.hh"
#include <iostream>
#include <algorithm>

namespace Amanzi {

void TimeStepManager::RegisterTimeEvent(double start, double period, double stop) {
  timeEvents_.push_back(TimeEvent(start, period, stop));
}

void TimeStepManager::RegisterTimeEvent(std::vector<double> times) {
  // make sure we only admit sorted arrays with unique entries
  std::vector<double> loc_times;
  loc_times = times;
  std::sort(loc_times.begin(), loc_times.end());
  loc_times.erase(std::unique(loc_times.begin(), loc_times.end(), near_equal), loc_times.end());  
  timeEvents_.push_back(TimeEvent(loc_times));
}

void TimeStepManager::RegisterTimeEvent(double time) {
  timeEvents_.push_back(TimeEvent(time));
}

  double TimeStepManager::TimeStep(double T, double dT) const {
  double next_T_all_events(1e99);
  // loop over all events to find the next event time
  for (std::list<TimeEvent>::const_iterator i = timeEvents_.begin(); i != timeEvents_.end(); ++i) {
    double next_T_this_event(1e99);
    // there are two possible types of events
    switch (i->Type()) {
      case TimeEvent::SPS: {
        if ( T < i->start_ ) {
          next_T_this_event = i->start_;
        } else if ( (i->stop_ == -1.0) || (T < i->stop_) ) {
          double n_periods = floor( (T - i->start_ )/i->period_ );
          double tmp = i->start_ + (n_periods+1.0)*i->period_;
	  if (tmp <= i->stop_) next_T_this_event = tmp;
        }
      }
      case TimeEvent::TIMES: {
        // we assume that i->times_ is ordered with unique elements       
        for (std::vector<double>::const_iterator j=i->times_.begin(); j!=i->times_.end(); ++j) {
          if (*j > T) {
            next_T_this_event = *j;
            break;
          }
        }
      }
    }
    next_T_all_events = std::min(next_T_all_events, next_T_this_event);
  }
  if (next_T_all_events == 1e99) return dT;
  double time_remaining(next_T_all_events - T);
  if (dT >= time_remaining) {
    return time_remaining;
  } else if ( dT > 0.75*time_remaining) {
    return 0.5*time_remaining;
  } else {
    return dT;
  } 
}

void TimeStepManager::print(std::ostream& os, double start, double end) const {
  // create a sorted array of the times between start and end and print it
  std::vector<double> print_times;
  for (std::list<TimeEvent>::const_iterator i = timeEvents_.begin(); i != timeEvents_.end(); ++i) {  
    switch (i->Type()) {
      case TimeEvent::SPS: {
        double time = i->start_;
        while (time <= end && time <= i->stop_) {
          print_times.push_back(time);
          time += i->period_;
        }
      }
      case TimeEvent::TIMES: {
        for (std::vector<double>::const_iterator j=i->times_.begin(); j!=i->times_.end(); ++j) {
          if (*j >= start && *j <= end) {
            print_times.push_back(*j);
          }
        }
      }
    }
  }
  std::sort(print_times.begin(), print_times.end());
  print_times.erase(std::unique(print_times.begin(), print_times.end(), near_equal), print_times.end());
  for (std::vector<double>::const_iterator i=print_times.begin(); i!=print_times.end(); ++i) {
    os << *i << " ";
  }
}
} // namespace Amanzi
