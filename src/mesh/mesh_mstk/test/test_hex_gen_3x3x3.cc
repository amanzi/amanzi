/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"


#include "AmanziMap.hh"
#include "AmanziComm.hh"

#include "MeshAudit.hh"

// Test generation of hex mesh in serial

TEST(MSTK_HEX_GEN_3x3x3)
{
  int i, j, k, err, nc, nf, nv;
  Amanzi::AmanziMesh::Set_ID faces[6], nodes[8];

  int NV = 64;
  int NF = 108;
  int NC = 27;

  auto comm = Amanzi::getDefaultComm();

  // Generate a mesh consisting of 3x3x3 elements

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK(
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, comm));

  nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NV, nv);

  nf = mesh->num_entities(Amanzi::AmanziMesh::FACE,
                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF, nf);

  nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC, nc);


  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> c2f("", 6);
  auto cell_map = mesh->cell_map(false);
  auto face_map = mesh->face_map(false);
  for (int c = cell_map->getMinLocalIndex(); c <= cell_map->getMaxLocalIndex();
       c++) {
    CHECK_EQUAL(cell_map->getGlobalElement(c),
                mesh->getGlobalElement(c, Amanzi::AmanziMesh::CELL));
    mesh->cell_get_faces(c, c2f);
    for (int j = 0; j < 6; j++) {
      int f = face_map->getLocalElement(c2f(j));
      CHECK(f == c2f(j));
    }
  }

  std::ofstream outfile("test/mstk_hex_gen_3x3x3.out");
  Amanzi::MeshAudit auditor(mesh, outfile);
  auditor.Verify();
}
