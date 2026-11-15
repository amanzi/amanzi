/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "UnitTest++.h"

#include "TimeStepManager.hh"
#include "VerboseObject.hh"


TEST(TIME_STEP_MANAGER)
{
  Amanzi::Utils::TimeStepManager TSM;

  std::vector<double> times;
  times.push_back(5.5);
  times.push_back(5.5); // create a duplicate entry
  times.push_back(1.0); // now times is not sorted by size
  TSM.RegisterTimeEvent(times);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(std::string("1 5.5 "), ss.str());
  }


  TSM.RegisterTimeEvent(1.2, 1.5, 100.0);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(std::string("1 1.2 2.7 4.2 5.5 5.7 7.2 8.7 "), ss.str());
  }

  TSM.RegisterTimeEvent(4.0, 0.5, 6.0);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(std::string("1 1.2 2.7 4 4.2 4.5 5 5.5 5.7 6 7.2 8.7 "), ss.str());
  }

  TSM.RegisterTimeEvent(5.0);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(std::string("1 1.2 2.7 4 4.2 4.5 5 5.5 5.7 6 7.2 8.7 "), ss.str());
  }

  // test the timestep lmiiter for a time that
  // is related to a time event defined by start-period-stop
  CHECK_CLOSE(0.1, TSM.TimeStep(4.1, .5), 1e-12);
  CHECK_CLOSE(0.05, TSM.TimeStep(1.0, 0.05), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(1.0, 0.1), 1e-12);
  CHECK_CLOSE(0.149, TSM.TimeStep(1.0, 0.149), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(1.0, 0.151), 1e-12);

  // test the timestep lmiiter for a time that
  // is related to a time event defined by and array of times
  CHECK_CLOSE(0.1, TSM.TimeStep(4.4, .5), 1e-12);
  CHECK_CLOSE(0.05, TSM.TimeStep(5.3, 0.05), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(5.3, 0.1), 1e-12);
  CHECK_CLOSE(0.149, TSM.TimeStep(5.3, 0.149), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(5.3, 0.151), 1e-12);

  // test the timestep limiter for a time that is larger than
  // any time specified in time events
  CHECK_CLOSE(100.0, TSM.TimeStep(200.0, 100.0), 1e-12);

  // test the timestep limiter for a time that is equal to
  // the largest time specified in time events
  CHECK_CLOSE(50.0, TSM.TimeStep(100.0, 50.0), 1e-12);
}


TEST(TIME_STEP_MANAGER_PRESCRIBED_MAXIMUM)
{
  Teuchos::ParameterList p;
  p.set("prescribed timesteps [s]", Teuchos::Array<double>{4., 6., 2.});
  p.set("prescribed timesteps are maximum", true);
  Amanzi::Utils::TimeStepManager manager(p);
  CHECK_CLOSE(1., manager.TimeStep(0., 1.), 1e-12);
  CHECK_CLOSE(3., manager.TimeStep(1., 20.), 1e-12);
  // A retry neither consumes an interval nor throws in maximum mode.
  CHECK_CLOSE(.25, manager.TimeStep(1., .25, true), 1e-12);
  CHECK_CLOSE(3., manager.TimeStep(1., 20.), 1e-12);
  CHECK_CLOSE(6., manager.TimeStep(4., 20.), 1e-12);
  CHECK_CLOSE(2., manager.TimeStep(10., 20.), 1e-12);
  // Once the schedule ends, ordinary adaptive/event control continues.
  CHECK_CLOSE(20., manager.TimeStep(12., 20.), 1e-12);
  Amanzi::Utils::TimeStepManager restarted(p);
  CHECK_CLOSE(3., restarted.TimeStep(7., 20.), 1e-12);
  restarted.RegisterTimeEvent(8.);
  CHECK_CLOSE(1., restarted.TimeStep(7., 20.), 1e-12);
}

TEST(TIME_STEP_MANAGER_PRESCRIBED_MAXIMUM_OFFSET)
{
  Teuchos::ParameterList p;
  p.set("prescribed timesteps [s]", Teuchos::Array<double>{4., 6.});
  p.set("prescribed timesteps are maximum", true);
  p.set("prescribed timesteps start time [s]", 5.);
  Amanzi::Utils::TimeStepManager manager(p);
  CHECK_CLOSE(4., manager.TimeStep(5., 20.), 1e-12);
  CHECK_CLOSE(6., manager.TimeStep(9., 20.), 1e-12);
}

TEST(TIME_STEP_MANAGER_STRICT_PRESCRIPTION_UNCHANGED)
{
  Teuchos::ParameterList p;
  p.set("prescribed timesteps [s]", Teuchos::Array<double>{4., 6.});
  Amanzi::Utils::TimeStepManager manager(p);
  CHECK_CLOSE(4., manager.TimeStep(0., 1.), 1e-12);
  CHECK_CLOSE(6., manager.TimeStep(4., 1.), 1e-12);
  bool threw = false;
  try { manager.TimeStep(10., 1., true); } catch (...) { threw = true; }
  CHECK(threw);
}


TEST(TIME_STEP_MANAGER_PRESCRIBED_MAXIMUM_ROUNDOFF)
{
  Teuchos::ParameterList p;
  p.set("prescribed timesteps [s]", Teuchos::Array<double>{4.});
  p.set("prescribed timesteps are maximum", true);
  p.set("prescribed timesteps start time [s]", 1.e7);
  Amanzi::Utils::TimeStepManager manager(p);
  CHECK_CLOSE(4., manager.TimeStep(1.e7, 4.-1.e-9), 1e-12);
}
