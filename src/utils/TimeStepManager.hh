/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Daniil Svyatskiy 
*/

#ifndef TIME_STEP_MANAGER__
#define TIME_STEP_MANAGER__

#include <cmath>
#include <list>
#include <ostream>
#include <vector>

#include "VerboseObject.hh"

namespace Amanzi {

bool inline near_equal(double x, double y)
{
  return (fabs(x - y) < 1e-12 * std::max(1.0, std::max(fabs(x), fabs(y))));
}

bool inline close(double x, double y, double eps)
{
  return (fabs(x - y) < eps * std::max(1.0, std::max(fabs(x), fabs(y))));
}

class TimeStepManager {
  struct TimeEvent {
    enum TimeType { SPS, TIMES };
    TimeEvent() {}
    TimeEvent(double start, double period, double stop, bool phys)
      : start_(start), period_(period), stop_(stop), type_(SPS), physical_(phys)
    {}
    TimeEvent(const std::vector<double>& times, bool phys)
      : times_(times), type_(TIMES), physical_(phys)
    {}
    explicit TimeEvent(double time, bool phys) : type_(TIMES), physical_(phys)
    {
      times_.push_back(time);
    }
    TimeType Type() const { return type_; }
    bool isPhysical() const { return physical_; }

    double start_, period_, stop_;
    std::vector<double> times_;
    TimeType type_;
    bool physical_; // physical means that this event effects time step due to change in physics
  };

 public:
  TimeStepManager();
  TimeStepManager(Teuchos::ParameterList& verb_list);
  TimeStepManager(Teuchos::RCP<VerboseObject> vo_cd);
  void RegisterTimeEvent(double start, double period, double stop, bool phys = true);
  void RegisterTimeEvent(std::vector<double> times, bool phys = true);
  void RegisterTimeEvent(double time, bool phys = true);
  double TimeStep(const double T, const double dT, bool after_failure = false);
  void print(std::ostream& os, double start, double end) const;

 private:
  std::list<TimeEvent> timeEvents_;
  double dt_stable_storage;

 protected:
  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace Amanzi

#endif // TIME_STEP_MANAGER__
