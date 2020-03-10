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

struct test_cv {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVectorSpace> x_space;
  Teuchos::RCP<CompositeVector> x;

  test_cv()
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
    x = x_space->Create();
  }
  ~test_cv() {}
};


SUITE(COMPOSITE_VECTOR)
{
  // test basic setup and construction
  TEST_FIXTURE(test_cv, CVConstruction)
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
        4, x->getMap()->ComponentMap("cell", false)->getNodeNumElements());
      CHECK_EQUAL(
        5, x->getMap()->ComponentMap("cell", true)->getNodeNumElements());
    } else {
      CHECK_EQUAL(
        8, x->getMap()->ComponentMap("cell", false)->getNodeNumElements());
      CHECK_EQUAL(
        8, x->getMap()->ComponentMap("cell", true)->getNodeNumElements());
    }

    // check zero initialization
    auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", false);
    CHECK_CLOSE(0.0, v1(0, 0), 0.00001);
  }

  TEST_FIXTURE(test_cv, CVSetGet)
  {
    // check putscalar
    x->putScalar(2.0);
    {
      // check on owned
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", false);
      CHECK_CLOSE(2.0, v1(0, 0), 0.00001);
      CHECK_CLOSE(2.0, v1(0, 1), 0.00001);

      auto v2 = x->ViewComponent<AmanziDefaultHost>("face", 0, false);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }


    {
      // check on ghosted
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", true);
      CHECK_CLOSE(2.0, v1(0, 0), 0.00001);
      CHECK_CLOSE(2.0, v1(0, 1), 0.00001);

      auto v2 = x->ViewComponent<AmanziDefaultHost>("face", 0, true);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }

    // PutScalar by component
    // Not yet implemented, but when it is, this should work.
    //
    // std::vector<double> vals(2);
    // vals[0] = 4.0; vals[1] = 5.0;
    // x->putScalar("cell", vals);
    // x->putScalar("face", 3.0);
    // {
    //   // check on owned
    //   auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", false);
    //   CHECK_CLOSE(4.0, v1(0,0), 0.00001);
    //   CHECK_CLOSE(5.0, v1(0,1), 0.00001);

    //   auto v2 = x->ViewComponent<AmanziDefaultHost>("face", false);
    //   CHECK_CLOSE(3.0, v2(0,0), 0.00001);
    // }
    // {
    //   // check on ghosted
    //   auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", true);
    //   CHECK_CLOSE(4.0, v1(0,0), 0.00001);
    //   CHECK_CLOSE(5.0, v1(0,1), 0.00001);

    //   auto v2 = x->ViewComponent<AmanziDefaultHost>("face", true);
    //   CHECK_CLOSE(3.0, v2(0,0), 0.00001);
    // }


    // check set by view
    {
      // set the value, then destroy the view
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      v1(0, 0) = 16.0;
    }
    {
      // check set by view on owned
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(16.0, v1(0), 0.00001);
    }
    {
      // check set by view on ghosted
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", 0, true);
      CHECK_CLOSE(16.0, v1(0), 0.00001);
    }
  }


  // test the vector's copy constructor
  TEST_FIXTURE(test_cv, CVCopy)
  {
    x->putScalar(2.0);

    CompositeVector y(*x);
    y.putScalar(4.0);

    {
      // check set by view on owned
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(2.0, v1(0), 0.00001);

      // check set by view on owned
      auto v2 = y.ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(4.0, v2(0), 0.00001);
    }
  }

  // test the vector's operator=
  TEST_FIXTURE(test_cv, CVOperatorEqual)
  {
    x->putScalar(2.0);

    CompositeVector y(*x);
    y.putScalar(0.0);

    // operator= and check vals
    y = *x;
    {
      // check set by view on owned
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(2.0, v1(0), 0.00001);

      // check set by view on owned
      auto v2 = y.ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }

    // ensure operator= did not copy pointers
    x->putScalar(4.0);
    {
      // check set by view on owned
      auto v1 = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(4.0, v1(0), 0.00001);

      // check set by view on owned
      auto v2 = y.ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }
  }

  // test the communication routines
  TEST_FIXTURE(test_cv, CVScatter)
  {
    int rank = comm->getRank();
    int size = comm->getSize();
    int ncells =
      mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    x->putScalar(rank + 1);
    x->ScatterMasterToGhosted("cell");

    if (size == 2) {
      auto x_c = x->ViewComponent<AmanziDefaultHost>("cell", 0, true);
      if (rank == 0) {
        CHECK_CLOSE(1.0, x_c(0), 0.00001);
        CHECK_CLOSE(2.0, x_c(4), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(2.0, x_c(0), 0.00001);
        CHECK_CLOSE(1.0, x_c(4), 0.00001);
      }
    } else {
      auto x_c = x->ViewComponent<AmanziDefaultHost>("cell", 0, true);
      CHECK_CLOSE(rank + 1, x_c(0), 0.00001);
    }
  }

  TEST_FIXTURE(test_cv, CVGather)
  {
    int rank = comm->getRank();
    int size = comm->getSize();
    int ncells =
      mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

    { // scope for x_c
      auto x_c = x->ViewComponent<AmanziDefaultHost>("cell", 0, true);
      for (int c = 0; c != ncells; ++c) { x_c(c) = rank + 1; }
    }
    x->GatherGhostedToMaster("cell");

    if (size == 2) {
      auto x_c = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      if (rank == 0) {
        CHECK_CLOSE(1.0, x_c(0), 0.00001);
        CHECK_CLOSE(3.0, x_c(3), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(3.0, x_c(0), 0.00001);
        CHECK_CLOSE(2.0, x_c(3), 0.00001);
      }
    } else {
      auto x_c = x->ViewComponent<AmanziDefaultHost>("cell", 0, false);
      CHECK_CLOSE(1, x_c(0), 0.00001);
      CHECK_CLOSE(1, x_c(7), 0.00001);
    }
  }
}
