#include <iostream>
#include "UnitTest++.h"

#include "time_step_manager.hh"

TEST(TIME_STEP_MANAGER) {

  Amanzi::TimeStepManager TSM;

  TSM.RegisterTimeEvent(1.2, 1.5, 100.0);

  CHECK_CLOSE(0.3, TSM.TimeStep(20.1, .5), 1e-12);

}
