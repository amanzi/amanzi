#include <cstdlib>
#include <string>

#include "UnitTest++.h"
#include "boost/algorithm/string.hpp"

#include "ChemistryUtilities.hh"

SUITE(amanzi_chemistry_unit_tests_ChemistryUtilities) {
  namespace acu = Amanzi::AmanziChemistry::utilities;

  TEST(TestUtilites_LnToLog) {
    double out = acu::ln_to_log(2.5);
    CHECK_CLOSE(out, 1.0857362047581294, 1.0e-10);
  }

  TEST(TestUtilities_LogToLn) {
    double out = acu::log_to_ln(2.5);
    CHECK_CLOSE(out, 5.7564627324851152, 1.0e-10);
  }
}
