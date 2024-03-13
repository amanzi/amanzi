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
#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"

int
main(int argc, char* argv[])
{
  int res = 0;
  {
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    Kokkos::initialize();
    res = UnitTest::RunAllTests();
    Kokkos::finalize();
  }
  return res;
}
