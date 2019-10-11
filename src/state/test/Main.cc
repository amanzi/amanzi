#include <mpi.h>

#include "Teuchos_GlobalMPISession.hpp"
#include <TestReporterStdout.h>
#include <UnitTest++.h>

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}
