#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

  MPI_Init(&argc,&argv);

  return UnitTest::RunAllTests ();

  MPI_Finalize();
}

