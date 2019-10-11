#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

#include "VerboseObject_objs.hh"

int
main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);

  int status = UnitTest::RunAllTests();

  Kokkos::finalize();
  MPI_Finalize();

  return status;
}
