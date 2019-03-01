#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "bilinear_form_registration.hh"
#include "state_evaluators_registration.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  return UnitTest::RunAllTests();
}

