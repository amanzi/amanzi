#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"
#include "wrm_flow_registration.hh"

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
