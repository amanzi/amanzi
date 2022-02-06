#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_transport_registration.hh"

int main( int argc, char *argv[] )
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  return UnitTest::RunAllTests();  
}

