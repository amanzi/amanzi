#include <UnitTest++.h>

#include "../dbc.hh"

#ifdef ENABLE_DBC

SUITE(DBC)
{
  TEST(Assert)
  {
    CHECK_THROW(AMANZI_ASSERT(1 == 2), DBC::Assertion);
  }
}

#endif
