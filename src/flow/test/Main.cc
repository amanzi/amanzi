#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"


#include "tests_flow_evaluator_reg.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  return UnitTest::RunAllTests();
}

