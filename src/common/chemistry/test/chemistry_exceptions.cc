#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include <UnitTest++.h>

#include "ChemistryException.hh"
#include "exceptions.hh"

SUITE(GeochemistryTests_ChemistryException) {
  namespace ac = Amanzi::AmanziChemistry;
  TEST(TestChemistryException_default_message) {
    try {
      ac::ChemistryException ce;
      Exceptions::amanzi_throw(ce);
    } catch (ac::ChemistryException& e) {
      CHECK_EQUAL("CHEMISTRY_ERROR: An unknown error has occured.", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryException_error_string) {
    CHECK_EQUAL("CHEMISTRY_ERROR: ", ac::kChemistryError);
  }  // end TEST()

  TEST(TestChemistryException_message) {
    try {
      ac::ChemistryException ce("Foo bar baz.");
      Exceptions::amanzi_throw(ce);
    } catch (ac::ChemistryException& e) {
      CHECK_EQUAL("CHEMISTRY_ERROR: Foo bar baz.", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryInvalidInput_message) {
    try {
      ac::ChemistryInvalidInput ce("Foo bar baz.");
      Exceptions::amanzi_throw(ce);
    } catch (ac::ChemistryException& e) {
      CHECK_EQUAL("CHEMISTRY_ERROR: Foo bar baz.", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryInvalidSolution_message) {
    try {
      ac::ChemistryInvalidSolution ce("Foo bar baz.");
      Exceptions::amanzi_throw(ce);
    } catch (ac::ChemistryException& e) {
      CHECK_EQUAL("CHEMISTRY_ERROR: Foo bar baz.", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryUnrecoverableError_message) {
    try {
      ac::ChemistryUnrecoverableError ce("Foo bar baz.");
      Exceptions::amanzi_throw(ce);
    } catch (ac::ChemistryException& e) {
      CHECK_EQUAL("CHEMISTRY_ERROR: Foo bar baz.", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryMaxIterationsReached_message) {
    try {
      ac::ChemistryMaxIterationsReached ce("Foo bar baz.");
      Exceptions::amanzi_throw(ce);
    } catch (ac::ChemistryException& e) {
      CHECK_EQUAL("CHEMISTRY_ERROR: Foo bar baz.", e.what());
    }
  }  // end TEST()
}  // end SUITE()
