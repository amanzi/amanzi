#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>
#include "Teuchos_GlobalMPISession.hpp"


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests ();
}

