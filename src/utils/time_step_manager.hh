#ifndef TIME_STEP_MANAGER__
#define TIME_STEP_MANAGER__

#include <vector>
#include <list>
#include <ostream>
#include <cmath>

namespace Amanzi {

bool inline near_equal (double x, double y) {
  return (fabs(x-y)<1e-12*std::max(1.0,std::max(fabs(x),fabs(y))));
}

class TimeStepManager {

  struct TimeEvent {
    enum TimeType { SPS, TIMES };
    TimeEvent() {}
    TimeEvent(double start, double period, double stop)
        : start_(start), period_(period), stop_(stop), type_(SPS) {}
    TimeEvent(const std::vector<double>& times)
        : times_(times), type_(TIMES) {}
    explicit TimeEvent(double time) : type_(TIMES) {
      times_.push_back(time);
    }
    TimeType Type() const {
      return type_;
    }

    double start_, period_, stop_;
    std::vector<double> times_;
    TimeType type_;
  };

 public:
  void RegisterTimeEvent(double start, double period, double stop);
  void RegisterTimeEvent(std::vector<double> times);
  void RegisterTimeEvent(double time);
  double TimeStep(const double T, const double dT) const;
  void print(std::ostream& os, double start, double end) const;

 private:
  std::list<TimeEvent> timeEvents_;

};
}

#endif // TIME_STEP_MANAGER__
