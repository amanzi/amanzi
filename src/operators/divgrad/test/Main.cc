#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>
#include "Teuchos_GlobalMPISession.hpp"

#include "test_mfd_surf.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  //  return UnitTest::RunAllTests ();
  test_mfd();
}

