/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include "../exceptions.hh"

using namespace Exceptions;

void
raise_amanzi_exception()
{
  Amanzi_exception exception;
  amanzi_throw(exception);
}

SUITE(Exceptions)
{
  TEST(Behavior)
  {
    // Check the default value of RAISE
    CHECK_EQUAL(exception_behavior(), RAISE);

    set_exception_behavior_abort();
    CHECK_EQUAL(exception_behavior(), ABORT);

    set_exception_behavior_raise();
    CHECK_EQUAL(exception_behavior(), RAISE);

    set_exception_behavior(ABORT);
    CHECK_EQUAL(exception_behavior(), ABORT);

    set_exception_behavior(RAISE);
    CHECK_EQUAL(exception_behavior(), RAISE);
  }

  TEST(Catch) { CHECK_THROW(raise_amanzi_exception(), Amanzi_exception); }
}
