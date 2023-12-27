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
}
