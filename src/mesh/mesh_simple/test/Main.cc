/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"

int
main(int argc, char* argv[])
{
  Kokkos::initialize();
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}
