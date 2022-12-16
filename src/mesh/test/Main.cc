#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

#include "VerboseObject_objs.hh"

#include "Kokkos_Core.hpp"

int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Kokkos::initialize();
  
  int status = UnitTest::RunAllTests ();
  
  Kokkos::finalize();
  MPI_Finalize();
  
  return status;
}

