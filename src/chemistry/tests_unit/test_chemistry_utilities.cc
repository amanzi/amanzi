/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>

#include <string>

#include <UnitTest++.h>

#include "chemistry_utilities.hh"

SUITE(amanzi_chemistry_unit_tests_ChemistryUtilities) {
  namespace acu = amanzi::chemistry::utilities;

  TEST(TestUtilities_CaseInsensitiveStringCompare) {
    std::string base_string("The quick brown fox jumps over the lazy dog.");
    std::string       test1("The Quick Brown Fox Jumps Over The Lazy Dog.");
    std::string       test3("This is a test. This is only a test. Failure");
    std::string       test2("This is a test.");
    CHECK(acu::CaseInsensitiveStringCompare(base_string, test1));
    CHECK(!acu::CaseInsensitiveStringCompare(base_string, test2));
    CHECK(!acu::CaseInsensitiveStringCompare(base_string, test3));
  }  // end TEST(CaseInsensitiveStringCompare)

  TEST(TestUtilities_CompareFabs) {
    double value1 = 3.14159;
    double value2 = -3.14159;
    double value3 = -2.71828;
    double value4 = 2.71827;
    double value5 = -2.71829;
    CHECK(!acu::CompareFabs(value1, value2));  // fabs are equal --> false
    CHECK(!acu::CompareFabs(value3, value4));  // 2.71828 > 2.71827 --> false
    CHECK(acu::CompareFabs(value4, value3));  // 2.71827 < 2.71828 --> true
    CHECK(!acu::CompareFabs(value5, value3));  // 2.71829 > 2.71828 --> false
    CHECK(acu::CompareFabs(value3, value5));  // 2.71828 < 2.71829 --> true
    CHECK(acu::CompareFabs(value4, value2));  // 2.71828 < 3.14159 --> true
  }  // end TEST(CompareFabs)

  TEST(TestUtilities_LowerCaseString) {
    std::string test_string("ThIS is A MiXeD CaSe STRING.");
    std::string expected("this is a mixed case string.");
    std::string received;
    acu::LowerCaseString(test_string, &received);
    CHECK_EQUAL(received, expected);
  }  // end TEST(LowerCaseString)

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
  }  // end TEST(RemoveLeadingAndTrailingWhitespace)

  TEST(TestUtilites_LnToLog) {
    double out = acu::ln_to_log(2.5);
    CHECK_CLOSE(out, 1.0857362047581294, 1.0e-10);
  }  // end TEST(LnToLog)

  TEST(TestUtilities_LogToLn) {
    double out = acu::log_to_ln(2.5);
    CHECK_CLOSE(out, 5.7564627324851152, 1.0e-10);
  }  // end TEST(LogToLn)

}  // end SUITE(ChemistryUtilities)
