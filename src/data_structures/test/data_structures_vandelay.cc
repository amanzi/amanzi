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

// #include "Mesh_MSTK.hh"
#include "AmanziComm.hh"
#include "MeshFactory.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_cv_vandelay {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<CompositeVectorSpace> x_space;
  Teuchos::RCP<CompositeVector> x;

  test_cv_vandelay()
  {
    comm = getDefaultComm();
    MeshFactory meshfactory(comm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    std::vector<Entity_kind> locations = { CELL, FACE };
    std::vector<std::string> names = { "cell", "face" };

    std::vector<int> num_dofs = { 2, 1 };

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);
    x = x_space->Create();
  }
  ~test_cv_vandelay() {}
};


SUITE(VANDELAY_COMPOSITE_VECTOR)
{
  // test the vector's putscalar
  TEST_FIXTURE(test_cv_vandelay, CVVandelay)
  {
    std::cout << "X has " << x->size() << " components" << std::endl;
    x->putScalar(2.0);

    {
      auto v_c = x->viewComponent<MirrorHost>("cell", false);
      CHECK_CLOSE(2., v_c(0, 0), 0.00001);
      CHECK_CLOSE(2., v_c(0, 1), 0.00001);

      auto v_f = x->viewComponent<MirrorHost>("face", 0, false);
      CHECK_CLOSE(2.0, v_f(0), 0.00001);
    }
    std::cout << " set componets are good" << std::endl;

    {
      auto v_bf = x->viewComponent<MirrorHost>("boundary_face", 0, false);
      CHECK_CLOSE(2.0, v_bf(0), 0.00001);
    }
  }

  // test owned vs ghosted
  TEST_FIXTURE(test_cv_vandelay, CVVandelayGhosted)
  {
    std::cout << "X has " << x->NumComponents() << " components" << std::endl;
    x->PutScalar(2.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);

    *x2->ViewComponent("boundary_face", false) = *x->ViewComponent("boundary_face", false);

    // test the scatter of boundary_faces
    x2->ScatterMasterToGhosted("boundary_face");
    int nbf_owned = x2->ViewComponent("boundary_face", false)->MyLength();
    int nbf_all = x2->ViewComponent("boundary_face", true)->MyLength();
    int size = comm->NumProc();
    if (size == 1) {
      CHECK_EQUAL(nbf_all, nbf_owned);
      CHECK_CLOSE((*x2)("boundary_face", 0, nbf_owned - 1), 2.0, 0.00001);
    } else {
      CHECK(nbf_owned < nbf_all);
      CHECK_CLOSE((*x2)("boundary_face", 0, nbf_all - 1), 2.0, 0.00001);
    }
  }
}
