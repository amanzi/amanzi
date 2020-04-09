#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"

int
main(int argc, char* argv[])
{
  Kokkos::initialize();
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int status = UnitTest::RunAllTests();

  Kokkos::finalize();
  return status;
}
