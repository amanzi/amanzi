/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "models_shallow_water_reg.hh"
#include "pks_shallow_water_reg.hh"
#include "VerboseObject_objs.hh"

#include "Kokkos_Core.hpp"

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  int status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}
