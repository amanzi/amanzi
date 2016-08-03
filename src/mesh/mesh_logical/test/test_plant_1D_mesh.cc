/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshLogicalFactory.hh"
#include "Geometry.hh"

#include "plant_1D_mesh.hh"

TEST(PLANT_1D_MESH_CONSTRUCTION)
{
  using namespace Amanzi;
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3));

  Teuchos::RCP<AmanziMesh::MeshLogical> m =
      Amanzi::Testing::plantMesh(&comm, gm, true);

  CHECK_EQUAL(1 + 10 + 10 + 10 + 50,
              m->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED));
  CHECK_EQUAL(1 + 10 + 10 + 10 + 10 + 51,
              m->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED));

  CHECK_EQUAL(1, m->get_set_size("leaf", AmanziMesh::CELL, AmanziMesh::OWNED));
  CHECK_EQUAL(10, m->get_set_size("stem", AmanziMesh::CELL, AmanziMesh::OWNED));
  CHECK_EQUAL(10, m->get_set_size("troot", AmanziMesh::CELL, AmanziMesh::OWNED));
  CHECK_EQUAL(10, m->get_set_size("aroot", AmanziMesh::CELL, AmanziMesh::OWNED));
  CHECK_EQUAL(50, m->get_set_size("soil", AmanziMesh::CELL, AmanziMesh::OWNED));

}
