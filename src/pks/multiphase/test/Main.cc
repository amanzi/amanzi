#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "eos_registration.hh"
#include "pks_multiphase_registration.hh"
#include "state_evaluators_registration.hh"
#include "wrmmp_registration.hh"

#include "VerboseObject_objs.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  return UnitTest::RunAllTests();
}

