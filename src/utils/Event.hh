/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#pragma once

namespace Amanzi {
namespace Utils {

//
// Either int or double
//
template<typename Scalar>
struct Event_EPS {
  static const int value = 0;
};

template<>
struct Event_EPS<double> {
  static const double value = 1.e-12;
};


//
// helper function to determine equality that may be either exact (int) or not (double)
//
template<typename Scalar>
bool inline isNearEqual(Scalar x, Scalar y, Scalar eps = Event_EPS<Scalar>::value)
{
  return (std::abs(x - y) <= eps * std::max((Scalar)1, std::max(std::abs(x), std::abs(y))));
};


//
// A class for storing sequences of events, as either int or doubles.
//
template<typename Scalar>
struct Event {
  Event() {}
  virtual Scalar getNext(Scalar time) const = 0;
  virtual bool contains(Scalar time) const = 0;
};


//
// A list of events
//
template<typename Scalar>
struct EventList : public Event<Scalar> {
  EventTimes(const std::vector<Scalar>& times)
    : times_(times)
  {
    // make sure we only admit sorted arrays with unique entries
    std::sort(times_.begin(), times_.end());
    times_.erase(std::unique(times_.begin(), times_.end(), isNearEqual), times_.end());
  }

  Scalar getNext(Scalar time) const override {
    for (auto j : times_) {
      if (j > time && !isNearEqual(j, time)) {
        return j;
      }
    }
    return (Scalar) -1;
  }

  bool contains(Scalar time) const override {
    return std::find_if(times_.begin(), times_.end(),
                        [=](Scalar j) { return isNearEqual(time, j) }) != times_.end();
  }

 private:
  std::vector<Scalar> times_;
};


//
// Events by (potentially open-ended) start, period, and stop
//
template<typename Scalar>
class EventSPS : public Event<Scalar> {
  EventSPS(Scalar start, Scalar period, Scalar stop)
    : start_(start),
      period_(period),
      stop_(stop)
  {
    if (period_ <= (Scalar)0) {
      Errors::Message msg("EventSPS provided invalid period.");
      Exceptions::amanzi_throw(msg);
    }
    if (stop_ > (Scalar)0 && stop_ <= start_) {
      Errors::Message msg("EventSPS provided invalid start/stop.");
      Exceptions::amanzi_throw(msg);
    }
    if (start_ < (Scalar)0) {
      Errors::Message msg("EventSPS provided invalid start.");
      Exceptions::amanzi_throw(msg);
    }
  }

  Scalar getNext(Scalar time) const override {
    if (time < start_ && !isNearEqual(time, start_)) {
      return start_;
    } else {
      int n_periods = static_cast<int>(std::floor((double)(time - start_) / period_));
      if (isNearEqual((double)(n_periods + 1), ((double)(time - start_) / period_))) n_periods += 1;
      Scalar tmp = start_ + (n_periods + 1) * period_;

      if (stop_ < 0) {
        return tmp;
      } else if (isNearEqual(tmp, stop_)) {
        return stop_;
      } else if (tmp > stop_) {
        return -1;
      } else {
        return tmp;
      }
    }
  }

  bool contains(Scalar time) const override {
    return isNearEqual((double)(time - start_) / period_, 0.);
  }

private:
  Scalar start_, period_, stop_;
};


} // namespace Utils
} // namespace Amanzi
