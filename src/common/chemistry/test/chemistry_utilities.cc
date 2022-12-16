/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <cstdlib>
#include <string>

#include "UnitTest++.h"

#include "ChemistryUtilities.hh"

SUITE(amanzi_chemistry_unit_tests_ChemistryUtilities)
{
  namespace acu = Amanzi::AmanziChemistry::utilities;

  TEST(TestUtilites_LnToLog)
  {
    double out = acu::ln_to_log(2.5);
    CHECK_CLOSE(out, 1.0857362047581294, 1.0e-10);
  }

  TEST(TestUtilities_LogToLn)
  {
    double out = acu::log_to_ln(2.5);
    CHECK_CLOSE(out, 5.7564627324851152, 1.0e-10);
  }
}
