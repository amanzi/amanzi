#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"
#include "wrm_flow_registration.hh"

int main(int argc, char *argv[])
{
  int res = 0;
  {
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    res = UnitTest::RunAllTests();
  }
  return res;
}

