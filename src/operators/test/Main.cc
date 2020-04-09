#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"

//#include "bilinear_form_registration.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc,argv);
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status; 
}

