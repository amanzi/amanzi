#include <cstdlib>
#include <string>

#include "UnitTest++.h"
#include "boost/algorithm/string.hpp"

#include "ChemistryUtilities.hh"

SUITE(amanzi_chemistry_unit_tests_ChemistryUtilities) {
  namespace acu = Amanzi::AmanziChemistry::utilities;

  TEST(TestUtilities_RemoveLeadingAndTrailingWhitespace) {
    std::string base("This is a test.");
    std::string leading(" \t This is a test.");
    std::string trailing("This is a test.\v \n");
    std::string lead_trail(" \f This is a test. \n\r ");
    acu::RemoveLeadingAndTrailingWhitespace(&leading);
    CHECK_EQUAL(base, leading);
    acu::RemoveLeadingAndTrailingWhitespace(&trailing);
    CHECK_EQUAL(base, trailing);
    acu::RemoveLeadingAndTrailingWhitespace(&lead_trail);
    CHECK_EQUAL(base, lead_trail);
  }

  TEST(TestUtilites_LnToLog) {
    double out = acu::ln_to_log(2.5);
    CHECK_CLOSE(out, 1.0857362047581294, 1.0e-10);
  }

  TEST(TestUtilities_LogToLn) {
    double out = acu::log_to_ln(2.5);
    CHECK_CLOSE(out, 5.7564627324851152, 1.0e-10);
  }
}
