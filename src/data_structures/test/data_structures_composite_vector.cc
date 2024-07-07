/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
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
    locations[0] = Entity_kind::CELL;
    locations[1] = Entity_kind::FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);
    x = x_space->Create();
  }
};


struct test_cv_bf {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVectorSpace> x_space;
  Teuchos::RCP<CompositeVector> x;

  test_cv_bf()
  {
    comm = getDefaultComm();
    MeshFactory meshfactory(comm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    std::vector<Entity_kind> locations(2);
    locations[0] = Entity_kind::CELL;
    locations[1] = Entity_kind::BOUNDARY_FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "boundary_face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);
    x = x_space->Create();
  }
};


SUITE(COMPOSITE_VECTOR)
{
  // test basic setup and construction
  TEST_FIXTURE(test_cv, CVConstruction)
  {
    CHECK_EQUAL(2, x->size());
    int size = comm->getSize();
    CHECK_EQUAL(2, x->getNumVectors("cell"));
    CHECK(x->getMap()->size() == x->size());
    CHECK(x->getMap()->isSameAs(*x->getMap()));
    CHECK(x->getMap()->getComponentMap("cell", false)->isSameAs(*mesh->getMap(Entity_kind::CELL, false)));
    CHECK(x->getMap()->getComponentMap("cell", true)->isSameAs(*mesh->getMap(Entity_kind::CELL, true)));

    if (size == 2) {
      CHECK_EQUAL(4, x->getMap()->getComponentMap("cell", false)->getLocalNumElements());
      CHECK_EQUAL(5, x->getMap()->getComponentMap("cell", true)->getLocalNumElements());
    } else {
      CHECK_EQUAL(8, x->getMap()->getComponentMap("cell", false)->getLocalNumElements());
      CHECK_EQUAL(8, x->getMap()->getComponentMap("cell", true)->getLocalNumElements());
    }

    // check zero initialization
    auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", false);
    CHECK_CLOSE(0.0, v1(0, 0), 0.00001);
  }

  TEST_FIXTURE(test_cv, CVSetGet)
  {
    // check putscalar
    x->putScalar(2.0);
    {
      // check on owned
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", false);
      CHECK_CLOSE(2.0, v1(0, 0), 0.00001);
      CHECK_CLOSE(2.0, v1(0, 1), 0.00001);

      auto v2 = x->viewComponent<MemSpace_kind::HOST>("face", 0, false);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }


    {
      // check on ghosted
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", true);
      CHECK_CLOSE(2.0, v1(0, 0), 0.00001);
      CHECK_CLOSE(2.0, v1(0, 1), 0.00001);

      auto v2 = x->viewComponent<MemSpace_kind::HOST>("face", 0, true);
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
    //   auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", false);
    //   CHECK_CLOSE(4.0, v1(0,0), 0.00001);
    //   CHECK_CLOSE(5.0, v1(0,1), 0.00001);

    //   auto v2 = x->viewComponent<MemSpace_kind::HOST>("face", false);
    //   CHECK_CLOSE(3.0, v2(0,0), 0.00001);
    // }
    // {
    //   // check on ghosted
    //   auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", true);
    //   CHECK_CLOSE(4.0, v1(0,0), 0.00001);
    //   CHECK_CLOSE(5.0, v1(0,1), 0.00001);

    //   auto v2 = x->viewComponent<MemSpace_kind::HOST>("face", true);
    //   CHECK_CLOSE(3.0, v2(0,0), 0.00001);
    // }


    // check set by view
    {
      // set the value, then destroy the view
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      v1(0) = 16.0;
    }
    {
      // check set by view on owned
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(16.0, v1(0), 0.00001);
    }
    {
      // check set by view on ghosted
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", 0, true);
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
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(2.0, v1(0), 0.00001);

      // check set by view on owned
      auto v2 = y.viewComponent<MemSpace_kind::HOST>("cell", 0, false);
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
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(2.0, v1(0), 0.00001);

      // check set by view on owned
      auto v2 = y.viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }

    // ensure operator= did not copy pointers
    x->putScalar(4.0);
    {
      // check set by view on owned
      auto v1 = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(4.0, v1(0), 0.00001);

      // check set by view on owned
      auto v2 = y.viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(2.0, v2(0), 0.00001);
    }
  }

  // test the communication routines
  TEST_FIXTURE(test_cv, CVScatter)
  {
    int rank = comm->getRank();
    int size = comm->getSize();
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    x->putScalar(rank + 1);
    x->scatterMasterToGhosted("cell");

    if (size == 2) {
      auto x_c = x->viewComponent<MemSpace_kind::HOST>("cell", 0, true);
      if (rank == 0) {
        CHECK_CLOSE(1.0, x_c(0), 0.00001);
        CHECK_CLOSE(2.0, x_c(4), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(2.0, x_c(0), 0.00001);
        CHECK_CLOSE(1.0, x_c(4), 0.00001);
      }
    } else {
      auto x_c = x->viewComponent<MemSpace_kind::HOST>("cell", 0, true);
      CHECK_CLOSE(rank + 1, x_c(0), 0.00001);
    }
  }

  TEST_FIXTURE(test_cv, CVGather)
  {
    int rank = comm->getRank();
    int size = comm->getSize();
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

    { // scope for x_c
      auto x_c = x->viewComponent<MemSpace_kind::HOST>("cell", 0, true);
      for (int c = 0; c != ncells; ++c) { x_c(c) = rank + 1; }
    }
    x->gatherGhostedToMaster("cell");

    if (size == 2) {
      auto x_c = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      if (rank == 0) {
        CHECK_CLOSE(1.0, x_c(0), 0.00001);
        CHECK_CLOSE(3.0, x_c(3), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(3.0, x_c(0), 0.00001);
        CHECK_CLOSE(2.0, x_c(3), 0.00001);
      }
    } else {
      auto x_c = x->viewComponent<MemSpace_kind::HOST>("cell", 0, false);
      CHECK_CLOSE(1, x_c(0), 0.00001);
      CHECK_CLOSE(1, x_c(7), 0.00001);
    }
  }

  TEST_FIXTURE(test_cv_bf, CVBoundaryFaces)
  {
    // this test just confirms that CVs use BOUNDARY_FACE maps correctly, a
    // long-running bug
    int size = comm->getSize();
    int rank = comm->getRank();

    int nbf_owned = x->viewComponent<MemSpace_kind::HOST>("boundary_face", false).extent(0);
    int nbf_all = x->viewComponent<MemSpace_kind::HOST>("boundary_face", true).extent(0);
    std::cout << "On rank " << rank << " of " << size << ", nbf_all = " << nbf_all
              << ", nbf_owned = " << nbf_owned << std::endl;

    if (size == 1) {
      CHECK_EQUAL(nbf_owned, nbf_all);
    } else {
      CHECK(nbf_owned < nbf_all);
    }
  }
}
