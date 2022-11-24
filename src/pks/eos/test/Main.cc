#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "eos_registration.hh"
#include "VerboseObject_objs.hh"

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
