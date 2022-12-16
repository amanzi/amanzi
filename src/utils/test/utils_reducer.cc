/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

#include "UnitTest++.h"
#include "StringReducer.hh"
#include "AmanziComm.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

SUITE(STRING_REDUCER)
{
  TEST(StringReducer_SAME)
  {
    auto comm = getDefaultComm();
    StringReducer<20> red(comm);

    // all the same
    std::vector<std::string> in{ "apples", "bananas", "canteloupes" };
    auto out = red.intersectAll(in);
    CHECK(out == in);
  }

  TEST(StringReducer_ALL_DIFFERENT)
  {
    auto comm = getDefaultComm();
    StringReducer<20> red(comm);

    int size = comm->NumProc();
    int rank = comm->MyPID();

    if (size > 1 && size < 7) {
      std::vector<std::string> all{ "apples", "bananas",      "canteloupes",
                                    "dates",  "elderberries", "figs" };
      std::vector<std::string> in{ all[rank] };
      auto out = red.intersectAll(in);
      CHECK_EQUAL(0, out.size());
    }
  }

  TEST(StringReducer_SOME_MISSING)
  {
    auto comm = getDefaultComm();
    StringReducer<20> red(comm);

    int size = comm->NumProc();
    int rank = comm->MyPID();

    if (size > 1) {
      // first and last are diff
      std::vector<std::string> in;
      if (rank == 1)
        in = { "apples", "bananas", "canteloupes", "dates", "elderberries" };
      else
        in = { "bananas", "dates" };

      auto out = red.intersectAll(in);
      CHECK_EQUAL(2, out.size());
      CHECK_EQUAL("bananas", out[0]);
      CHECK_EQUAL("dates", out[1]);
    }
  }

  TEST(StringReducer_DIFFERENT_MISSING)
  {
    auto comm = getDefaultComm();
    StringReducer<20> red(comm);

    int size = comm->NumProc();
    int rank = comm->MyPID();

    if (size > 1) {
      // first and last are diff
      std::vector<std::string> in;
      if (rank == 1)
        in = { "apples", "bananas", "dates" };
      else
        in = { "apples", "canteloupes", "dates" };

      auto out = red.intersectAll(in);
      CHECK_EQUAL(2, out.size());
      CHECK_EQUAL("apples", out[0]);
      CHECK_EQUAL("dates", out[1]);
    }
  }


  TEST(StringReducer_SOME_EMPTY)
  {
    auto comm = getDefaultComm();
    StringReducer<20> red(comm);

    // int size = comm->NumProc();
    int rank = comm->MyPID();

    // first and last are diff
    std::vector<std::string> in;
    if (rank != 0) in = { "apples", "bananas", "dates" };

    auto out = red.intersectAll(in);
    CHECK_EQUAL(0, out.size());
  }

} // SUITE
