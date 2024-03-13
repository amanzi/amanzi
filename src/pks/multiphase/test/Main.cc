/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "eos_reg.hh"
#include "evaluators_multiphase_reg.hh"
#include "models_energy_reg.hh"
#include "models_multiphase_reg.hh"
#include "pks_multiphase_reg.hh"
#include "state_evaluators_registration.hh"

#include "VerboseObject_objs.hh"

// a custom model for testing
#include "WRMmp_Custom.hh"

namespace Amanzi {
namespace Multiphase {

Utils::RegisteredFactory<WRMmp, WRMmp_Custom> WRMmp_Custom::reg_("Custom");

} // namespace Multiphase
} // namespace Amanzi

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  int status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}
