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
#include "CompositeVector.hh"
#include "TreeVector.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_tv {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVectorSpace> x_vec_space;
  Teuchos::RCP<CompositeVector> x_vec;
  Teuchos::RCP<TreeVector> x;
  Teuchos::RCP<TreeVector> x2;

  test_tv()
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

    x_vec_space = Teuchos::rcp(new CompositeVectorSpace());
    x_vec_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);
    x_vec = Teuchos::rcp(new CompositeVector(*x_vec_space));
    x = Teuchos::rcp(new TreeVector());
    x->SetData(x_vec);

    x2 = Teuchos::rcp(new TreeVector());
    x2->PushBack(x);
    x2->PushBack(x);
  }
  ~test_tv() {}
};


SUITE(TREE_VECTOR)
{
  // test the vector's putscalar
  TEST_FIXTURE(test_tv, TVPutScalar)
  {
    x->PutScalar(2.0);
    CHECK_CLOSE((*x->Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("face", 0, 0), 2.0, 0.00001);

    x2->PutScalar(3.0);
    CHECK_CLOSE((*x->Data())("cell", 0, 0),
                3.0,
                0.00001); // x2 created via PushBack(x), ensure this stores a pointer, not copy
    CHECK_CLOSE((*x->Data())("cell", 1, 0), 3.0, 0.00001);
    CHECK_CLOSE((*x->Data())("face", 0, 0), 3.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 0, 0), 3.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 1, 0), 3.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("face", 0, 0), 3.0, 0.00001);
  }

  // test the vector's copy constructor
  TEST_FIXTURE(test_tv, TVCopy)
  {
    x->PutScalar(2.0);

    TreeVector y(*x);
    CHECK(y.Map().SameAs(x->Map()));
    y.PutScalar(4.0);
    CHECK_CLOSE((*x->Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y.Data())("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*y.Data())("cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE((*y.Data())("face", 0, 0), 4.0, 0.00001);

    TreeVector z(x->Map());
    CHECK(z.Map().SameAs(x->Map()));
    z.PutScalar(4.0);
    CHECK_CLOSE((*x->Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*z.Data())("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*z.Data())("cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE((*z.Data())("face", 0, 0), 4.0, 0.00001);

    x2->PutScalar(2.0);
    TreeVector y2(*x2);
    CHECK(y2.Map().SameAs(x2->Map()));
    y2.PutScalar(5.0);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("face", 0, 0), 5.0, 0.00001);

    TreeVector z2(x2->Map());
    CHECK(z2.Map().SameAs(x2->Map()));
    z2.PutScalar(6.0);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*z2.SubVector(0)->Data())("cell", 0, 0), 6.0, 0.00001);
    CHECK_CLOSE((*z2.SubVector(0)->Data())("cell", 1, 0), 6.0, 0.00001);
    CHECK_CLOSE((*z2.SubVector(0)->Data())("face", 0, 0), 6.0, 0.00001);
  }


  // test the vector's operator=
  TEST_FIXTURE(test_tv, TVOperatorEqual)
  {
    x->PutScalar(2.0);

    TreeVector y(*x);
    y.PutScalar(0.0);

    // operator= and check vals
    y = *x;
    CHECK_CLOSE((*x->Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x->Data())("face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y.Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y.Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y.Data())("face", 0, 0), 2.0, 0.00001);

    // ensure operator= did not copy pointers
    x->PutScalar(4.0);
    CHECK_CLOSE((*x->Data())("cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*x->Data())("cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE((*x->Data())("face", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE((*y.Data())("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y.Data())("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*y.Data())("face", 0, 0), 2.0, 0.00001);

    x2->PutScalar(5.0);
    TreeVector y2(*x2);
    y2.PutScalar(0.0);

    // operator= and check vals
    y2 = *x2;

    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("face", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("face", 0, 0), 5.0, 0.00001);

    // ensure operator= did not copy pointers
    x2->PutScalar(6.0);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 0, 0), 6.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("cell", 1, 0), 6.0, 0.00001);
    CHECK_CLOSE((*x2->SubVector(0)->Data())("face", 0, 0), 6.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE((*y2.SubVector(0)->Data())("face", 0, 0), 5.0, 0.00001);
  }
}
