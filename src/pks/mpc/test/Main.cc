#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "state_evaluators_registration.hh"
#include "constitutive_relations_eos_registration.hh"
#include "flow_relations_registration.hh"
#include "flow_permafrost_registration.hh"
#include "flow_richards_registration.hh"
#include "energy_constitutive_relations_internal_energy_registration.hh"
#include "energy_constitutive_relations_source_terms_registration.hh"
#include "energy_constitutive_relations_thermal_conductivity_registration.hh"
#include "energy_two_phase_registration.hh"
#include "energy_three_phase_registration.hh"
#include "mpc_registration.hh"


//#include "test_mpc_subsurface.cc"
int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return UnitTest::RunAllTests();
//  run();
}

