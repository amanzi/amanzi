/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Tests for state as a container of meshes
*/

// TPLs
#include "UnitTest++.h"

// State
#include "MeshFactory.hh"
#include "State.hh"

TEST(STATE_CREATION) {
  using namespace Amanzi;

  auto comm = Amanzi::getDefaultComm();
  AmanziMesh::MeshFactory fac(comm);
  auto mesh = fac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  State s;
  s.RegisterDomainMesh(mesh);

  s.RegisterMesh("other_domain_mesh", mesh, true);

  CHECK(s.HasMesh("domain"));
  CHECK(s.HasMesh("other_domain_mesh"));

  CHECK(!s.IsDeformableMesh("domain"));
  CHECK(s.IsDeformableMesh("other_domain_mesh"));
}
