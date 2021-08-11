#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>
#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"

// a fake WRM model for testing surface-subsurface coupling
#include "WRM_fake.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRM_fake> WRM_fake::factory_("fake");

}  // namespace Flow
}  // namespace Amanzi

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  
  return UnitTest::RunAllTests ();
}

