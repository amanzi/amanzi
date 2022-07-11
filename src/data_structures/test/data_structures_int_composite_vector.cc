/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_int_cv {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVectorSpace> x_space;
  Teuchos::RCP<CompositeVector_<int>> x;

  test_int_cv()
  {
    comm = getDefaultComm();
    MeshFactory meshfactory(comm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 8, 1, 1);

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(
      names, locations, num_dofs);
    x = x_space->Create<int>();
  }
  ~test_int_cv() {}
};


SUITE(COMPOSITE_VECTOR_INT)
{
  // test basic setup and construction
  TEST_FIXTURE(test_int_cv, CVConstruction)
  {
    CHECK_EQUAL(2, x->size());
    int size = comm->getSize();
    if (size == 1) CHECK_EQUAL(8, x->getLocalLength("cell"));
    CHECK_EQUAL(2, x->getNumVectors("cell"));
    CHECK(x->getMap()->size() == x->size());
    CHECK(x->getMap()->SameAs(*x->getMap()));
    CHECK(x->getMap()
            ->ComponentMap("cell", false)
            ->isSameAs(*mesh->map(CELL, false)));
    CHECK(x->getMap()
            ->ComponentMap("cell", true)
            ->isSameAs(*mesh->map(CELL, true)));

    if (size == 2) {
      CHECK_EQUAL(
        4, x->getMap()->ComponentMap("cell", false)->getLocalNumElements());
      CHECK_EQUAL(
        5, x->getMap()->ComponentMap("cell", true)->getLocalNumElements());
    } else {
      CHECK_EQUAL(
        8, x->getMap()->ComponentMap("cell", false)->getLocalNumElements());
      CHECK_EQUAL(
        8, x->getMap()->ComponentMap("cell", true)->getLocalNumElements());
    }
  }

  TEST_FIXTURE(test_int_cv, CVSetGet)
  {
    // check putscalar
    x->putScalar(2);
    {
      // check on owned
      auto v1 = x->ViewComponent<MirrorHost>("cell", false);
      CHECK_EQUAL(2, v1(0, 0));
      CHECK_EQUAL(2, v1(0, 1));

      auto v2 = x->ViewComponent<MirrorHost>("face", 0, false);
      CHECK_EQUAL(2, v2(0));
    }


    {
      // check on ghosted
      auto v1 = x->ViewComponent<MirrorHost>("cell", true);
      CHECK_EQUAL(2, v1(0, 0));
      CHECK_EQUAL(2, v1(0, 1));

      auto v2 = x->ViewComponent<MirrorHost>("face", 0, true);
      CHECK_EQUAL(2, v2(0));
    }

    // PutScalar by component
    // Not yet implemented, but when it is, this should work.
    //
    // std::vector<double> vals(2);
    // vals[0] = 4; vals[1] = 5;
    // x->putScalar("cell", vals);
    // x->putScalar("face", 3);
    // {
    //   // check on owned
    //   auto v1 = x->ViewComponent<MirrorHost>("cell", false);
    //   CHECK_EQUAL(4, v1(0,0));
    //   CHECK_EQUAL(5, v1(0,1));

    //   auto v2 = x->ViewComponent<MirrorHost>("face", false);
    //   CHECK_EQUAL(3, v2(0,0));
    // }
    // {
    //   // check on ghosted
    //   auto v1 = x->ViewComponent<MirrorHost>("cell", true);
    //   CHECK_EQUAL(4, v1(0,0));
    //   CHECK_EQUAL(5, v1(0,1));

    //   auto v2 = x->ViewComponent<MirrorHost>("face", true);
    //   CHECK_EQUAL(3, v2(0,0));
    // }


    // check set by view
    {
      // set the value, then destroy the view
      auto v1 = x->ViewComponent<MirrorHost>("cell", 0, false);
      v1(0) = 16;
    }
    {
      // check set by view on owned
      auto v1 = x->ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(16, v1(0));
    }
    {
      // check set by view on ghosted
      auto v1 = x->ViewComponent<MirrorHost>("cell", 0, true);
      CHECK_EQUAL(16, v1(0));
    }
  }


  // test the vector's copy constructor
  TEST_FIXTURE(test_int_cv, CVCopy)
  {
    x->putScalar(2);

    CompositeVector_<int> y(*x);
    y.putScalar(4);

    {
      // check set by view on owned
      auto v1 = x->ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(2, v1(0));

      // check set by view on owned
      auto v2 = y.ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(4, v2(0));
    }
  }

  // test the vector's operator=
  TEST_FIXTURE(test_int_cv, CVOperatorEqual)
  {
    x->putScalar(2);

    CompositeVector_<int> y(*x);
    y.putScalar(0);

    // operator= and check vals
    y = *x;
    {
      // check set by view on owned
      auto v1 = x->ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(2, v1(0));

      // check set by view on owned
      auto v2 = y.ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(2, v2(0));
    }

    // ensure operator= did not copy pointers
    x->putScalar(4);
    {
      // check set by view on owned
      auto v1 = x->ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(4, v1(0));

      // check set by view on owned
      auto v2 = y.ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(2, v2(0));
    }
  }

  // test the communication routines
  TEST_FIXTURE(test_int_cv, CVScatter)
  {
    int rank = comm->getRank();
    int size = comm->getSize();
    int ncells =
      mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    x->putScalar(rank + 1);
    x->ScatterMasterToGhosted("cell");

    if (size == 2) {
      auto x_c = x->ViewComponent<MirrorHost>("cell", 0, true);
      if (rank == 0) {
        CHECK_EQUAL(1, x_c(0));
        CHECK_EQUAL(2, x_c(4));
      } else if (rank == 1) {
        CHECK_EQUAL(2, x_c(0));
        CHECK_EQUAL(1, x_c(4));
      }
    } else {
      auto x_c = x->ViewComponent<MirrorHost>("cell", 0, true);
      CHECK_EQUAL(rank + 1, x_c(0));
    }
  }

  TEST_FIXTURE(test_int_cv, CVGather)
  {
    int rank = comm->getRank();
    int size = comm->getSize();
    int ncells =
      mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

    { // scope for x_c
      auto x_c = x->ViewComponent<MirrorHost>("cell", 0, true);
      for (int c = 0; c != ncells; ++c) { x_c(c) = rank + 1; }
    }
    x->GatherGhostedToMaster("cell");

    if (size == 2) {
      auto x_c = x->ViewComponent<MirrorHost>("cell", 0, false);
      if (rank == 0) {
        CHECK_EQUAL(1, x_c(0));
        CHECK_EQUAL(3, x_c(3));
      } else if (rank == 1) {
        CHECK_EQUAL(3, x_c(0));
        CHECK_EQUAL(2, x_c(3));
      }
    } else {
      auto x_c = x->ViewComponent<MirrorHost>("cell", 0, false);
      CHECK_EQUAL(1, x_c(0));
      CHECK_EQUAL(1, x_c(7));
    }
  }
}
