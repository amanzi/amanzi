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

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "composite_vector.hh"
#include "tree_vector.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_tv {
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVector> x_vec;
  Teuchos::RCP<TreeVector> x;

  test_tv() {
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

    x_vec = Teuchos::rcp(new CompositeVector(mesh, names, locations, num_dofs, true));
    x_vec->CreateData();
    x = Teuchos::rcp(new TreeVector("x"));
    x->set_data(x_vec);
  }
  ~test_tv() { delete comm; }
};


SUITE(TREE_VECTOR) {
  // test the vector's putscalar
  TEST_FIXTURE(test_tv, TVPutScalar) {
    x->PutScalar(2.0);
    CHECK_CLOSE((*x->data())("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x->data())("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*x->data())("face",0,0), 2.0, 0.00001);
  }

  // test the vector's copy constructor
  TEST_FIXTURE(test_tv, TVCopy) {
    x->PutScalar(2.0);

    TreeVector y(*x);
    y.PutScalar(4.0);
    CHECK_CLOSE((*x->data())("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x->data())("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*x->data())("face",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*y.data())("cell",0,0), 4.0, 0.00001);
    CHECK_CLOSE((*y.data())("cell",1,0), 4.0, 0.00001);
    CHECK_CLOSE((*y.data())("face",0,0), 4.0, 0.00001);
  }

  // test the vector's operator=
  TEST_FIXTURE(test_tv, TVOperatorEqual) {
    x->PutScalar(2.0);

    TreeVector y(*x);
    y.PutScalar(0.0);

    // operator= and check vals
    y = *x;
    CHECK_CLOSE((*x->data())("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*x->data())("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*x->data())("face",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*y.data())("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*y.data())("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*y.data())("face",0,0), 2.0, 0.00001);

    // ensure operator= did not copy pointers
    x->PutScalar(4.0);
    CHECK_CLOSE((*x->data())("cell",0,0), 4.0, 0.00001);
    CHECK_CLOSE((*x->data())("cell",1,0), 4.0, 0.00001);
    CHECK_CLOSE((*x->data())("face",0,0), 4.0, 0.00001);
    CHECK_CLOSE((*y.data())("cell",0,0), 2.0, 0.00001);
    CHECK_CLOSE((*y.data())("cell",1,0), 2.0, 0.00001);
    CHECK_CLOSE((*y.data())("face",0,0), 2.0, 0.00001);
  }
}

