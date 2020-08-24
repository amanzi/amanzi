#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  return UnitTest::RunAllTests();
}

