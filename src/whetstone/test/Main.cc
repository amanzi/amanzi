#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"
#include "bilinear_form_registration.hh"


int main( int argc, char *argv[] )
{
  Teuchos::GlobalMPISession  mpiSession( &argc, &argv );

  return UnitTest::RunAllTests();  
}

