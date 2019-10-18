/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <mpi.h>

#include "Teuchos_GlobalMPISession.hpp"
#include <TestReporterStdout.h>
#include <UnitTest++.h>

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"

#include "test/factory_models_reg.hh"

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}
