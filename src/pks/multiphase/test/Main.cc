/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>

#include "Teuchos_GlobalMPISession.hpp"

#include "eos_registration.hh"
#include "energy_iem_registration.hh"
#include "energy_tcm_registration.hh"
#include "multiphase_evaluators_registration.hh"
#include "pks_multiphase_registration.hh"
#include "state_evaluators_registration.hh"
#include "wrmmp_registration.hh"

#include "VerboseObject_objs.hh"

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  return UnitTest::RunAllTests();
}
