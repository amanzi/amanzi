/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <ThrowingTestReporter.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();

  UnitTest::TestReporterStdout base_reporter;
  UnitTest::ThrowingTestReporter reporter(&base_reporter);
  UnitTest::TestRunner runner(reporter);
  int status = runner.RunTestsIf(UnitTest::Test::GetTestList(), NULL, UnitTest::True(), 0);
  Kokkos::finalize();
  return status;
}
