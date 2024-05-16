#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "numerical_flux_registration.hh"
#include "pks_shallow_water_registration.hh"
#include "VerboseObject_objs.hh"

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  return UnitTest::RunAllTests();
}
