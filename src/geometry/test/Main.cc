/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

#include "VerboseObject_objs.hh"

int
main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int status = UnitTest::RunAllTests();

  MPI_Finalize();

  return status;
}
