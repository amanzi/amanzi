/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Ethan Coon

   Unit tests for the composite vector.
   ------------------------------------------------------------------------- */

#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "composite_vector.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_cv {
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVector> x;
  Teuchos::RCP<CompositeVector> x2;

  test_cv() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x = Teuchos::rcp(new CompositeVector(mesh, names, locations, num_dofs, true));
    //    x2 = Teuchos::rcp(new CompositeVector(mesh, CELL, 1, true));
  }
  ~test_cv() { delete comm; }
};


SUITE(COMPOSITE_VECTOR) {
  // test the vector's putscalar
  TEST_FIXTURE(test_cv, CVPutScalar) {
    x->CreateData();
    x->PutScalar(2.0);
    CHECK_CLOSE((*x)("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 2.0, 0.00001);

    std::vector<double> vals(2);
    vals[0] = 4.0; vals[1] = 5.0;
    x->PutScalar("cell", vals);
    x->PutScalar("face", 3.0);
    CHECK_CLOSE((*x)("cell",0,0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("cell",1,0), 5.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 3.0, 0.00001);
  }

  // test the vector's set value/operator()
  TEST_FIXTURE(test_cv, CVSetValue) {
    x->CreateData();
    (*x)("cell",0,0) = 2.0;
    (*x)("cell",1,0) = 3.0;
    (*x)("face",0,0) = 4.0;

    CHECK_CLOSE((*x)("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell",1,0), 3.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 4.0, 0.00001);
  }

  // test the vector's view mechanisms
  TEST_FIXTURE(test_cv, CVViews) {
    x->CreateData();
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
  TEST_FIXTURE(test_cv, CVCopy) {
    x->CreateData();
    x->PutScalar(2.0);

    CompositeVector y(*x);
    y.PutScalar(4.0);
    CHECK_CLOSE((*x)("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell",0,0), 4.0, 0.00001);
    CHECK_CLOSE(y("cell",1,0), 4.0, 0.00001);
    CHECK_CLOSE(y("face",0,0), 4.0, 0.00001);
  }

  // test the vector's operator=
  TEST_FIXTURE(test_cv, CVOperatorEqual) {
    x->CreateData();
    x->PutScalar(2.0);

    CompositeVector y(*x);
    y.PutScalar(0.0);

    // operator= and check vals
    y = *x;
    CHECK_CLOSE((*x)("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE(y("face",0,0), 2.0, 0.00001);

    // ensure operator= did not copy pointers
    x->PutScalar(4.0);
    CHECK_CLOSE((*x)("cell",0,0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("cell",1,0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 4.0, 0.00001);
    CHECK_CLOSE(y("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE(y("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE(y("face",0,0), 2.0, 0.00001);
  }

  // test the vector's view of local vals
  TEST_FIXTURE(test_cv, CVView) {
    x->CreateData();
    x->PutScalar(2.0);

    // get a non-ghosted view
    Teuchos::RCP<Epetra_MultiVector> cells = x->ViewComponent("cell", false);
    CHECK_CLOSE((*cells)[0][0], 2.0, 0.00001);

    // alter the view, communicate, and check the ghosted view
    cells->PutScalar(4.0);
    x->ScatterMasterToGhosted();
    CHECK_CLOSE((*x)("cell",0,0), 4.0, 0.00001);
    CHECK_CLOSE((*x)("face",0,0), 2.0, 0.00001);
  }

  // DEPRECATED FEATURE
  // // test the vector's conversion to a packed SuperVector
  // TEST_FIXTURE(test_cv, CVSuperVector) {
  //   x->CreateData();
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
  TEST_FIXTURE(test_cv, CVScatter) {
    x->CreateData();
    int rank = comm->MyPID();
    int size = comm->NumProc();
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int c=0; c!=ncells; ++c) {
      (*x)("cell",0,c) = rank+1.0;
    }
    x->ScatterMasterToGhosted("cell");

    if (size == 2) {
      if (rank == 0) {
        CHECK_CLOSE(1.0, (*x)("cell",0), 0.00001);
        CHECK_CLOSE(2.0, (*x)("cell",5), 0.00001);
      } else if (rank == 1) {
        CHECK_CLOSE(2.0, (*x)("cell",0), 0.00001);
        CHECK_CLOSE(1.0, (*x)("cell",5), 0.00001);
      }
    } else {
      std::cout << "Test CVScatter requires 2 procs" << std::endl;
    }

    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    for (int f=0; f!=nfaces; ++f) {
      (*x)("face",0,f) = rank+1.0;
    }
    x->ScatterMasterToGhosted("face");
  }

  TEST_FIXTURE(test_cv, CVGather) {
    x->CreateData();
    int rank = comm->MyPID();
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    for (int c=0; c!=ncells; ++c) {
      (*x)("cell",0,c) = rank+1.0;
    }
    x->GatherGhostedToMaster("cell");

    int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
    for (int f=0; f!=nfaces; ++f) {
      (*x)("face",0,f) = rank+1.0;
    }
    x->GatherGhostedToMaster("face");
  }
}

