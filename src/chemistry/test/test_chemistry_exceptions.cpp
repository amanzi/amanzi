/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#include "UnitTest++.h"

#include "ChemistryException.hpp"

SUITE(GeochemistryTests_ChemistryException)
{
  TEST(TestChemistryException_default_message)
  {
    try {
      ChemistryException ce;
      throw ce;
    }
    catch (ChemistryException& e) {
      CHECK_EQUAL("ERROR", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryException_message)
  {
    try {
      ChemistryException ce("Foo bar baz.");
      throw ce;
    }
    catch (ChemistryException& e) {
      CHECK_EQUAL("Foo bar baz.", e.what());
    }
  }  // end TEST()

  TEST(TestChemistryException_status)
  {
    try {
      ChemistryException ce("Foo bar baz.", ChemistryException::kUnrecoverableError);
      throw ce;
    }
    catch (ChemistryException& e) {
      CHECK_EQUAL(ChemistryException::kUnrecoverableError, e.error_status());
    }
  }  // end TEST()


}  // end SUITE()


