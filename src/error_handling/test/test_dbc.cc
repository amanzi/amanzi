#include <UnitTest++.h>

#include "../dbc.hh"

#ifdef ENABLE_DBC

SUITE (DBC)
{

    TEST (Assert)
    {
        CHECK_THROW (ASSERT (1==2), DBC::Assertion);
    }

}

#endif
