/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include "UnitTest++.h"

#include "TimeStepManager.hh"
#include "VerboseObject.hh"


TEST(TIME_STEP_MANAGER)
{
  Amanzi::TimeStepManager TSM;

  std::vector<double> times;
  times.push_back(5.5);
  times.push_back(5.5); // create a duplicate entry
  times.push_back(1.0); // now times is not sorted by size
  TSM.RegisterTimeEvent(times);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(ss.str(), std::string("1 5.5 "));
  }


  TSM.RegisterTimeEvent(1.2, 1.5, 100.0);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(ss.str(), std::string("1 1.2 2.7 4.2 5.5 5.7 7.2 8.7 "));
  }

  TSM.RegisterTimeEvent(4.0, 0.5, 6.0);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(ss.str(),
                std::string("1 1.2 2.7 4 4.2 4.5 5 5.5 5.7 6 7.2 8.7 "));
  }

  TSM.RegisterTimeEvent(5.0);

  {
    std::stringstream ss;
    TSM.print(ss, 0.0, 10.0);
    CHECK_EQUAL(ss.str(),
                std::string("1 1.2 2.7 4 4.2 4.5 5 5.5 5.7 6 7.2 8.7 "));
  }

  // test the time step lmiiter for a time that
  // is related to a time event defined by start-period-stop
  CHECK_CLOSE(0.1, TSM.TimeStep(4.1, .5), 1e-12);
  CHECK_CLOSE(0.05, TSM.TimeStep(1.0, 0.05), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(1.0, 0.1), 1e-12);
  CHECK_CLOSE(0.149, TSM.TimeStep(1.0, 0.149), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(1.0, 0.151), 1e-12);

  // test the time step lmiiter for a time that
  // is related to a time event defined by and array of times
  CHECK_CLOSE(0.1, TSM.TimeStep(4.4, .5), 1e-12);
  CHECK_CLOSE(0.05, TSM.TimeStep(5.3, 0.05), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(5.3, 0.1), 1e-12);
  CHECK_CLOSE(0.149, TSM.TimeStep(5.3, 0.149), 1e-12);
  CHECK_CLOSE(0.1, TSM.TimeStep(5.3, 0.151), 1e-12);

  // test the time step limiter for a time that is larger than
  // any time specified in time events
  CHECK_CLOSE(100.0, TSM.TimeStep(200.0, 100.0), 1e-12);

  // test the time step limiter for a time that is equal to
  // the largest time specified in time events
  CHECK_CLOSE(50.0, TSM.TimeStep(100.0, 50.0), 1e-12);
}
