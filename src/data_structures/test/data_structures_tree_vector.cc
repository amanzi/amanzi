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
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
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
    x_vec_space->SetMesh(mesh)->SetGhosted()->SetComponents(
      names, locations, num_dofs);
    x_vec = x_vec_space->Create();

    auto x_space = Teuchos::rcp(new TreeVectorSpace(comm));
    x_space->SetData(x_vec->getMap());

    x = Teuchos::rcp(new TreeVector(x_space));
    x->SetData(x_vec);

    auto x2_space = Teuchos::rcp(new TreeVectorSpace(comm));
    x2_space->PushBack(x_space);
    x2_space->PushBack(x_space);
    x2 = Teuchos::rcp(new TreeVector(x2_space));
    x2->SetSubVector(0, x);
    x2->SetSubVector(1, x);
  }
  ~test_tv() {}
};


// NOTE: this is backwards notation, with dof_num first, as was the case in
// Epetra.
double
get_value(const CompositeVector& cv, const std::string& cname, int dof_num,
          int lid)
{
  auto vec = cv.ViewComponent<DefaultHost>(cname, true);
  return vec(lid, dof_num);
}


SUITE(TREE_VECTOR)
{
  TEST_FIXTURE(test_tv, TVDefaultZero)
  {
    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0), 0.0, 1.e-10);
  }


  // test the vector's putscalar
  TEST_FIXTURE(test_tv, TVPutScalar)
  {
    x->putScalar(2.0);

    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "face", 0, 0), 2.0, 0.00001);

    x2->putScalar(3.0);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0),
                3.0,
                0.00001); // x2 created via PushBack(x), ensure this stores a
                          // pointer, not copy
    CHECK_CLOSE(get_value(*x->Data(), "cell", 1, 0), 3.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "face", 0, 0), 3.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 0, 0), 3.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 1, 0), 3.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "face", 0, 0), 3.0, 0.00001);
  }

  // test the vector's copy constructor
  TEST_FIXTURE(test_tv, TVCopy)
  {
    x->putScalar(2.0);

    TreeVector y(*x);
    CHECK(y.getMap()->SameAs(*x->getMap()));
    y.putScalar(4.0);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "face", 0, 0), 4.0, 0.00001);

    TreeVector z(x->getMap());
    CHECK(z.getMap()->SameAs(*x->getMap()));
    z.putScalar(4.0);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*z.Data(), "cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*z.Data(), "cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*z.Data(), "face", 0, 0), 4.0, 0.00001);

    x2->putScalar(2.0);
    TreeVector y2(*x2);
    CHECK(y2.getMap()->SameAs(*x2->getMap()));
    y2.putScalar(5.0);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "face", 0, 0), 5.0, 0.00001);

    TreeVector z2(x2->getMap());
    CHECK(z2.getMap()->SameAs(*x2->getMap()));
    z2.putScalar(6.0);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(
      get_value(*z2.SubVector(0)->Data(), "cell", 0, 0), 6.0, 0.00001);
    CHECK_CLOSE(
      get_value(*z2.SubVector(0)->Data(), "cell", 1, 0), 6.0, 0.00001);
    CHECK_CLOSE(
      get_value(*z2.SubVector(0)->Data(), "face", 0, 0), 6.0, 0.00001);
  }


  // test the vector's operator=
  TEST_FIXTURE(test_tv, TVOperatorEqual)
  {
    x->putScalar(2.0);

    TreeVector y(*x);
    y.putScalar(0.0);

    // operator= and check vals
    y = *x;
    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "face", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "face", 0, 0), 2.0, 0.00001);

    // ensure operator= did not copy pointers
    x->putScalar(4.0);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "cell", 1, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*x->Data(), "face", 0, 0), 4.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE(get_value(*y.Data(), "face", 0, 0), 2.0, 0.00001);

    x2->putScalar(5.0);
    TreeVector y2(*x2);
    y2.putScalar(0.0);

    // operator= and check vals
    y2 = *x2;

    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "face", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "face", 0, 0), 5.0, 0.00001);

    // ensure operator= did not copy pointers
    x2->putScalar(6.0);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 0, 0), 6.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "cell", 1, 0), 6.0, 0.00001);
    CHECK_CLOSE(
      get_value(*x2->SubVector(0)->Data(), "face", 0, 0), 6.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "cell", 0, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "cell", 1, 0), 5.0, 0.00001);
    CHECK_CLOSE(
      get_value(*y2.SubVector(0)->Data(), "face", 0, 0), 5.0, 0.00001);
  }
}
