/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

#include "UnitTest++.h"
#include "StringReducer.hh"
#include "Reductions.hh"
#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"

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

    int size = comm->getSize();
    int rank = comm->getRank();

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

    int size = comm->getSize();
    int rank = comm->getRank();

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

    int size = comm->getSize();
    int rank = comm->getRank();

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

    // int size = comm->getSize();
    int rank = comm->getRank();

    // first and last are diff
    std::vector<std::string> in;
    if (rank != 0) in = { "apples", "bananas", "dates" };

    auto out = red.intersectAll(in);
    CHECK_EQUAL(0, out.size());
  }

} // SUITE


TEST(VECTOR_REDUCTIONS)
{
  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(3 * comm->getSize(), 3, 0, comm));
  Vector_type vec(map);

  int my_rank = comm->getRank();
  {
    auto vv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
    vv(0, 0) = my_rank + 1;
    vv(1, 0) = my_rank + 2;
    vv(2, 0) = my_rank;
  }

  auto min_loc = Amanzi::Reductions::reduceAllMinLoc(vec);
  auto max_loc = Amanzi::Reductions::reduceAllMaxLoc(vec);
  CHECK_EQUAL(0.0, min_loc.val);
  CHECK_EQUAL(comm->getSize() - 1 + 2, max_loc.val);
  if (comm->getSize() == 1) {
    // rank 0, lid 2
    CHECK_EQUAL(2, min_loc.loc);
    // rank 0, lid 1
    CHECK_EQUAL(1, max_loc.loc);
  } else if (comm->getSize() == 2) {
    // rank 0, lid 2
    CHECK_EQUAL(2, min_loc.loc);
    // rank 1, lid 1
    CHECK_EQUAL(4, max_loc.loc);
  } else if (comm->getSize() == 5) {
    // rank 0, lid 2
    CHECK_EQUAL(2, min_loc.loc);
    // rank 4, lid 1
    CHECK_EQUAL(13, max_loc.loc);
  } else {
    AMANZI_ASSERT(false);
  }
}
