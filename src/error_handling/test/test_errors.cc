/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>
#include <iostream>

#include "../errors.hh"

SUITE(Messages)
{
  TEST(Add_data)
  {
    Errors::Message m("a");
    m.add_data("b");
    m << "c";

    CHECK_EQUAL("abc", m.what());
  }

  TEST(Copy)
  {
    const char* test_string = "abcdefghijklmnopqrstuxwxyz";

    Errors::Message m(test_string);
    Errors::Message m2(m);

    CHECK_EQUAL(test_string, m2.what());
  }


  TEST(Throw)
  {
    try {
      Errors::Message message("Something went wrong: ");
      message << "FOO!=BAR";
      Exceptions::amanzi_throw(message);
    } catch (const Errors::Message& m) {
      CHECK_EQUAL("Something went wrong: FOO!=BAR", m.what());
    }
  }

  TEST(Multilevel)
  {
    try {
      try {
        Errors::Message message("FAILED.");
        CHECK_EQUAL("FAILED.", message.what());
        Exceptions::amanzi_throw(message);
      } catch (Errors::Message& m) {
        m << "CAUGHT.";
        CHECK_EQUAL("FAILED.CAUGHT.", m.what());
        throw m;
      }
    } catch (Exceptions::Amanzi_exception& e) {
      CHECK_EQUAL("FAILED.CAUGHT.", e.what());
    }
  }
}
