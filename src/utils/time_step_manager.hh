#ifndef TIME_STEP_MANAGER__
#define TIME_STEP_MANAGER__

#include <vector>
#include <list>

namespace Amanzi {

class TimeStepManager {

  struct TimeEvent {
    enum TimeType { SPS, TIMES };
    TimeEvent() {}
    TimeEvent(double start, double period, double stop)
        : start_(start), period_(period), stop_(stop), type_(SPS) {}
    TimeEvent(const std::vector<double>& cycles) 
        : cycles_(cycles), type_(TIMES) {}
    TimeType Type() const { return type_; }

    double start_, period_, stop_;
    std::vector<double> cycles_;
    TimeType type_;
  };  

 public:
  void RegisterTimeEvent(double start, double period, double stop);
  void RegisterTimeEvent(const std::vector<double>& cycles);
  double TimeStep(const double T, const double dT);

 private:
  std::list<TimeEvent> timeEvents_;

};
}

#endif // TIME_STEP_MANAGER__
