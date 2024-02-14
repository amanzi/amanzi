/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"

#include "Kokkos_Core.hpp"

int
main(int argc, char* argv[])
{
  int res;
  {
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    Kokkos::initialize();
    {
      res = UnitTest::RunAllTests();
    }
    Kokkos::finalize();
  }
  return res;
}
