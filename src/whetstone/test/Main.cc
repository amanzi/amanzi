#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"

#include "bilinear_form_registration.hh"


int
main(int argc, char* argv[])
{
  Kokkos::initialize();
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int status = UnitTest::RunAllTests();

  Kokkos::finalize();
  return status;
}
