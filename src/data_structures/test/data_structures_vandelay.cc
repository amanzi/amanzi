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

// #include "Mesh_MSTK.hh"
#include "MeshFactory.hh"
#include "CompositeVector.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_cv_vandelay {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVectorSpace> x_space;
  Teuchos::RCP<CompositeVectorSpace> x2_space;
  Teuchos::RCP<CompositeVector> x;
  Teuchos::RCP<CompositeVector> x2;

  test_cv_vandelay()
  {
    comm = getDefaultComm();
    Preference pref;
    pref.clear();
    pref.push_back(Framework::MSTK);

    MeshFactory meshfactory(comm);
    meshfactory.set_preference(pref);

    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    // mesh = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2, comm, NULL));

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<Entity_kind> locations_v(2);
    locations_v[0] = CELL;
    locations_v[1] = BOUNDARY_FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<std::string> names_v(2);
    names_v[0] = "cell";
    names_v[1] = "boundary_face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(names, locations, num_dofs);

    x2_space = Teuchos::rcp(new CompositeVectorSpace());
    x2_space->SetMesh(mesh)->SetGhosted()->SetComponents(names_v, locations_v, num_dofs);

    x = Teuchos::rcp(new CompositeVector(*x_space));
    x2 = Teuchos::rcp(new CompositeVector(*x2_space));
  }
  ~test_cv_vandelay() {}
};


SUITE(VANDELAY_COMPOSITE_VECTOR)
{
  // test the vector's putscalar
  TEST_FIXTURE(test_cv_vandelay, CVVandelay)
  {
    std::cout << "X has " << x->NumComponents() << " components" << std::endl;
    x->PutScalar(2.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);

    *x2->ViewComponent("boundary_face", false) =
      *std::as_const(*x).ViewComponent("boundary_face", false);
    CHECK_CLOSE((*x2)("boundary_face", 0, 0), 2.0, 0.00001);
  }

  // test owned vs ghosted
  TEST_FIXTURE(test_cv_vandelay, CVVandelayGhosted)
  {
    std::cout << "X has " << x->NumComponents() << " components" << std::endl;
    x->PutScalar(2.0);
    CHECK_CLOSE((*x)("cell", 0, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("cell", 1, 0), 2.0, 0.00001);
    CHECK_CLOSE((*x)("face", 0, 0), 2.0, 0.00001);

    *x2->ViewComponent("boundary_face", false) =
      *std::as_const(*x).ViewComponent("boundary_face", false);

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
