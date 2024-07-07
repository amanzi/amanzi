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
  Teuchos::RCP<CompositeVectorSpace> x2_space;
  Teuchos::RCP<CompositeVector> x2;

  test_cv_vandelay()
  {
    comm = getDefaultComm();
    MeshFactory meshfactory(comm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    std::vector<Entity_kind> locations = { Entity_kind::CELL, Entity_kind::FACE };
    std::vector<Entity_kind> locations_v = { Entity_kind::CELL, Entity_kind::BOUNDARY_FACE };
    std::vector<std::string> names = { "cell", "face" };
    std::vector<std::string> names_v = { "cell", "boundary_face" };


    std::vector<int> num_dofs = { 2, 1 };

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);
    x = x_space->Create();

    x2_space = Teuchos::rcp(new CompositeVectorSpace());
    x2_space->SetMesh(mesh)->SetGhosted()->SetComponents(names_v, locations_v, num_dofs);
    x2 = x2_space->Create();
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
      auto v_c = x->viewComponent<MemSpace_kind::HOST>("cell", false);
      CHECK_CLOSE(2., v_c(0, 0), 0.00001);
      CHECK_CLOSE(2., v_c(0, 1), 0.00001);

      auto v_f = x->viewComponent<MemSpace_kind::HOST>("face", 0, false);
      CHECK_CLOSE(2.0, v_f(0), 0.00001);
    }
    std::cout << " set componets are good" << std::endl;

    {
      auto v_bf = x->viewComponent<MemSpace_kind::HOST>("boundary_face", 0, false);
      CHECK_CLOSE(2.0, v_bf(0), 0.00001);
    }
  }

  // test owned vs ghosted
  TEST_FIXTURE(test_cv_vandelay, CVVandelayGhosted)
  {
    std::cout << "X has " << x->size() << " components" << std::endl;
    x->putScalar(2.0);
    {
      auto v_c = x->viewComponent<MemSpace_kind::HOST>("cell", false);
      CHECK_CLOSE(2.0, v_c(0, 0), 0.00001);
      CHECK_CLOSE(2.0, v_c(0, 1), 0.00001);

      auto v_f = x->viewComponent<MemSpace_kind::HOST>("face", false);
      CHECK_CLOSE(2.0, v_f(0, 0), 0.00001);
    }

    Kokkos::deep_copy(x2->viewComponent("boundary_face", false),
                      x->viewComponent("boundary_face", false));

    // test the scatter of boundary_faces
    x2->scatterMasterToGhosted("boundary_face");
    int nbf_owned = x2->getComponent("boundary_face", false)->getLocalLength();
    int nbf_all = x2->getComponent("boundary_face", true)->getLocalLength();
    int size = comm->getSize();
    if (size == 1) {
      CHECK_EQUAL(nbf_all, nbf_owned);
      auto v_bf = x2->viewComponent<MemSpace_kind::HOST>("boundary_face", false);
      CHECK_CLOSE(v_bf(nbf_owned - 1, 0), 2.0, 0.00001);
    } else {
      CHECK(nbf_owned < nbf_all);
      auto v_bf = x2->viewComponent<MemSpace_kind::HOST>("boundary_face", false);
      CHECK_CLOSE(2.0, v_bf(nbf_all - 1, 0), 0.00001);
    }
  }
}
