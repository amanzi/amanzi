/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "UnitTest++.h"

#include "state_evaluators_registration.hh"
#include "VerboseObject_objs.hh"

#include "evaluators_flow_reg.hh"
#include "models_flow_reg.hh"

// a fake model for testing
#include "WRM_fake.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRM_fake> WRM_fake::factory_("fake");

} // namespace Flow
} // namespace Amanzi

int
main(int argc, char* argv[])
{
  int res = 0;
  {
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    res = UnitTest::RunAllTests();
  }
  return res;
}
