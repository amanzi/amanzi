#ifndef TIME_STEP_MANAGER__
#define TIME_STEP_MANAGER__

#include <vector>
#include <list>
#include <ostream>

namespace Amanzi {

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
  void RegisterTimeEvent(const std::vector<double>& cycles);
  void RegisterTimeEvent(double time);
  double TimeStep(const double T, const double dT) const;
  void print(std::ostream& os, double start, double end) const;

 private:
  std::list<TimeEvent> timeEvents_;

};
}

#endif // TIME_STEP_MANAGER__
