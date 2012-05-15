#include <iostream>
#include "UnitTest++.h"

#include "time_step_manager.hh"

TEST(TIME_STEP_MANAGER) {

  Amanzi::TimeStepManager TSM;

  TSM.RegisterTimeEvent(1.2, 1.5, 100.0);

  CHECK_EQUAL(0.75, TSM.TimeStep(20.1, .5));

}
