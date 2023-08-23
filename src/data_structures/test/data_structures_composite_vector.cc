/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   Amanzi

   Unit tests for the composite vector.
------------------------------------------------------------------------- */

#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "CompositeVector.hh"

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
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

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
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);
    x = Teuchos::rcp(new CompositeVector(*x_space));
  }
  ~test_cv() {}
};


SUITE(COMPOSITE_VECTOR)
{
  // test basic setup and construction
  TEST_FIXTURE(test_cv, CVConstruction)
  {
    CompositeVectorSpace space(x->Map());
    CHECK(x->Map().SameAs(space));
    CHECK_EQUAL(2, x->NumComponents());
    int size = comm->NumProc();
    if (size == 1) CHECK_EQUAL(8, x->size("cell"));
    CHECK_EQUAL(2, x->NumVectors("cell"));
    CHECK(x->Map().NumComponents() == x->NumComponents());
    CHECK(x->Map().SameAs(x->Map()));
  }

  // test the vector's putscalar
  TEST_FIXTURE(test_cv, CVPutScalar)
  {
    x->PutScalar(2.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);

    std::vector<double> vals(2);
    vals[0] = 4.0;
    vals[1] = 5.0;
    x->PutScalar("cell", vals);
    x->PutScalar("face", 3.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 3.0, 0.00001);
  }

  // test the vector's set value/operator()
  TEST_FIXTURE(test_cv, CVSetValue)
  {
    Epetra_MultiVector& x_c = *x->ViewComponent("cell", false);
    x_c[0][0] = 2.;
    x_c[1][0] = 3.;
    Epetra_MultiVector& x_f = *x->ViewComponent("face", false);
    x_f[0][0] = 4;

    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 3.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 4.0, 0.00001);
  }

  // test the vector's view mechanisms
  TEST_FIXTURE(test_cv, CVViews)
  {
    x->PutScalar(2.0);

    // check the values got set correctly and we can access them from views
    CHECK_CLOSE((*x->ViewComponent("cell", true))[0][0], 2.0, 0.00001);
    CHECK_CLOSE((*x->ViewComponent("cell", false))[0][0], 2.0, 0.00001);
    CHECK_CLOSE((*x->ViewComponent("cell"))[0][0], 2.0, 0.00001);

    // change the non-ghosted view and check the ghosted view
    x->ViewComponent("cell", false)->PutScalar(3.0);
    CHECK_CLOSE((*x->ViewComponent("cell", true))[0][0], 3.0, 0.00001);
    CHECK_CLOSE((*x->ViewComponent("cell", false))[0][0], 3.0, 0.00001);
    CHECK_CLOSE((*x->ViewComponent("cell"))[0][0], 3.0, 0.00001);
  }

  // test the vector's copy constructor
  TEST_FIXTURE(test_cv, CVCopy)
  {
    x->PutScalar(2.0);

    CompositeVector y(*x);
    y.PutScalar(4.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE(y("cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE(y("face", 0, 0), 4.0, 0.00001);
  }

  // test the vector's operator=
  TEST_FIXTURE(test_cv, CVOperatorEqual)
  {
    x->PutScalar(2.0);

    CompositeVector y(*x);
    y.PutScalar(0.0);

    // operator= and check vals
    y = *x;
    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(y("face", 0, 0), 2.0, 0.00001);

    // ensure operator= did not copy pointers
    x->PutScalar(4.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE(y("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(y("face", 0, 0), 2.0, 0.00001);
  }

  // test the vector's view of local vals
  TEST_FIXTURE(test_cv, CVView)
  {
    x->PutScalar(2.0);

    // get a non-ghosted view
    Teuchos::RCP<Epetra_MultiVector> cells = x->ViewComponent("cell", false);
    CHECK_CLOSE((*cells)[0][0], 2.0, 0.00001);

    // alter the view, communicate, and check the ghosted view
    cells->PutScalar(4.0);
    x->ScatterMasterToGhosted();
    CHECK_CLOSE((*x)("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);
  }

  // DEPRECATED FEATURE
  // // test the vector's conversion to a packed SuperVector
  // TEST_FIXTURE(test_cv, CVSuperVector) {
  //   x->ViewComponent("cell", true)->PutScalar(2.0);
  //   x->ViewComponent("face", true)->PutScalar(3.0);

  //   // create a supervector and check the values
  //   Teuchos::RCP<Epetra_Vector> y = x->CreateSuperVector();
  //   x->CopyToSuperVector(y);
  //   CHECK_CLOSE((*y)[0], 2.0, 0.00001);
  //   CHECK_CLOSE((*y)[43], 3.0, 0.00001);

  //   // copy construct the supervector and get new values back
  //   Epetra_Vector z(*y);
  //   z.PutScalar(5.0);
  //   x->DataFromSuperVector(z);
  //   CHECK_CLOSE((*x)("cell",0,0), 5.0, 0.00001);
  //   CHECK_CLOSE((*x)("face",0,0), 5.0, 0.00001);
  // }

  // test the communication routines
  TEST_FIXTURE(test_cv, CVScatter)
  {
    int rank = comm->MyPID();
    int size = comm->NumProc();
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

    { // scope for x_c
      Epetra_MultiVector& x_c = *x->ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) { x_c[0][c] = rank + 1.0; }
      x->ScatterMasterToGhosted("cell");
    }

    if (size == 2) {
      if (rank == 0) {
        CHECK_CLOSE(1.0, (*x)("cell", 0), 0.00001);
        CHECK_CLOSE(2.0, (*x)("cell", 5), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(2.0, (*x)("cell", 0), 0.00001);
        CHECK_CLOSE(1.0, (*x)("cell", 5), 0.00001);
      }
    } else {
      std::cout << "Test CVScatter requires 2 procs" << std::endl;
    }

    { // scope for x_f
      int nfaces =
        mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      Epetra_MultiVector& x_f = *x->ViewComponent("face", false);
      for (int f = 0; f != nfaces; ++f) { x_f[0][f] = rank + 1.0; }
      x->ScatterMasterToGhosted("face");
    }
  }

  TEST_FIXTURE(test_cv, CVGather)
  {
    int rank = comm->MyPID();
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
    Epetra_MultiVector& x_c = *x->ViewComponent("cell", true);
    for (int c = 0; c != ncells; ++c) { x_c[0][c] = rank + 1; }
    x->GatherGhostedToMaster("cell");

    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);
    Epetra_MultiVector& x_f = *x->ViewComponent("face", true);
    for (int f = 0; f != nfaces; ++f) { x_f[0][f] = rank + 1; }
    x->GatherGhostedToMaster("face");
  }

  TEST_FIXTURE(test_cv, CVManageCommSmoke)
  {
    // Ensures that Communication happens after a change.
    int rank = comm->MyPID();
    int size = comm->NumProc();
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
    Epetra_MultiVector& x_c = *x->ViewComponent("cell", false);
    for (int c = 0; c != ncells; ++c) { x_c[0][c] = rank + 1.0; }
    x->ScatterMasterToGhosted("cell");

    x->PutScalar(100.);
    x->ScatterMasterToGhosted("cell");

    if (size == 2) {
      if (rank == 0) {
        CHECK_CLOSE(100.0, (*x)("cell", 0), 0.00001);
        CHECK_CLOSE(100.0, (*x)("cell", 5), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(100.0, (*x)("cell", 0), 0.00001);
        CHECK_CLOSE(100.0, (*x)("cell", 5), 0.00001);
      }
    } else {
      std::cout << "Test CVScatter requires 2 procs" << std::endl;
    }
  }

  TEST_FIXTURE(test_cv, CVManageCommNoComm)
  {
    // Ensures that Communication DOESN"T happen if it isn't needed
    // Does not work
    int rank = comm->MyPID();
    int size = comm->NumProc();
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

    Epetra_MultiVector& x_c = *x->ViewComponent("cell", false);
    for (int c = 0; c != ncells; ++c) { x_c[0][c] = rank + 1.0; }
    x->ScatterMasterToGhosted("cell");
    if (size == 2) {
      if (rank == 0) {
        CHECK_CLOSE(1.0, (*x)("cell", 0), 0.00001);
        CHECK_CLOSE(2.0, (*x)("cell", 5), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(2.0, (*x)("cell", 0), 0.00001);
        CHECK_CLOSE(1.0, (*x)("cell", 5), 0.00001);
      }
    } else {
      std::cout << "Test CVScatter requires 2 procs" << std::endl;
    }

    // by changing the values using the same vector ref, we can cheat the
    // system
    for (int c = 0; c != ncells; ++c) { x_c[0][c] = rank + 100; }
    x->ScatterMasterToGhosted("cell"); // this call should NOT communicate
    if (size == 2) {
      if (rank == 0) {
        CHECK_CLOSE(100.0, (*x)("cell", 0), 0.00001);
        // CHECK_CLOSE(2.0, (*x)("cell",5), 0.00001);
        CHECK_CLOSE(101.0, (*x)("cell", 5), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(101.0, (*x)("cell", 0), 0.00001);
        // CHECK_CLOSE(1.0, (*x)("cell",5), 0.00001);
        CHECK_CLOSE(100.0, (*x)("cell", 5), 0.00001);
      }
    } else {
      std::cout << "Test CVScatter requires 2 procs" << std::endl;
    }
  }
}
