#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

  MPI_Init(&argc,&argv);

  int status = UnitTest::RunAllTests ();
  
  MPI_Finalize();
  
  return status;
}

