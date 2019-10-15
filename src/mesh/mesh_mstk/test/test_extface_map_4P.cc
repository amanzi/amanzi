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
#include "MeshAudit.hh"

#include "AmanziMap.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"


TEST(MSTK_EXTFACE_MAP_4P)
{
  int i, j, k, err, nc, nf, nv;
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<int> facedirs(6);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  auto comm = Amanzi::getDefaultComm();
  int rank = comm->getRank();
  int size = comm->getSize();
  CHECK_EQUAL(4, size);

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(
    new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_3x3x3_sets.exo", comm));

  auto face_map = mesh->face_map(false);
  auto extface_map = mesh->exterior_face_map(false);

  auto all_to_extface_importer = mesh->exterior_face_importer();

  for (int f = extface_map->getMinLocalIndex();
       f <= extface_map->getMaxLocalIndex();
       f++) {
    int gid = extface_map->getGlobalElement(f);
    int f2 = face_map->getLocalElement(gid); // f2 is local face id in face_map

    CHECK_EQUAL(face_map->getGlobalElement(f2), gid);

    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> fcells;
    mesh->face_get_cells(f2, Amanzi::AmanziMesh::Parallel_type::OWNED, fcells);
    CHECK_EQUAL(1, fcells.extent(0));
  }

  Amanzi::Vector_type allvec(face_map);
  Amanzi::Vector_type bdryvec(extface_map);

  // Insert the GlobalID of each face offsetted by 3 into the allvec
  {
    auto allvec_v = allvec.getDataNonConst();
    for (int f = face_map->getMinLocalIndex(); f < face_map->getMaxLocalIndex();
         f++)
      allvec_v[f] = face_map->getGlobalElement(f) + 3;
  }

  bdryvec.doImport(
    allvec, *all_to_extface_importer, Tpetra::CombineMode::INSERT);

  // Check if the importer got the right values from allvec into bdryvec
  // by checking if the values in the bdryvec minus the offset correspond
  // to the correct global IDs.
  bdryvec.sync_host();
  {
    auto bdryvec_v = bdryvec.getData();
    for (int f = extface_map->getMinLocalIndex();
         f < extface_map->getMaxLocalIndex();
         f++)
      CHECK_EQUAL(extface_map->getGlobalElement(f), bdryvec_v[f] - 3);
  }

  // Check if ghostmap contains only boundary faces
  auto extface_map_wghost = mesh->exterior_face_map(true);
  int nowned_bnd = extface_map->getNodeNumElements();
  int nnotowned_bnd = extface_map_wghost->getNodeNumElements() - nowned_bnd;

  std::vector<int> gl_id(nnotowned_bnd), pr_id(nnotowned_bnd),
    lc_id(nnotowned_bnd);

  for (int f = 0; f < nnotowned_bnd; f++) {
    gl_id[f] = extface_map_wghost->getGlobalElement(f + nowned_bnd);
  }

  auto gl_id_av = Teuchos::arrayView(gl_id.data(), nnotowned_bnd);
  auto pr_id_av = Teuchos::arrayView(pr_id.data(), nnotowned_bnd);
  auto lc_id_av = Teuchos::arrayView(lc_id.data(), nnotowned_bnd);
  extface_map->getRemoteIndexList(gl_id_av, pr_id_av, lc_id_av);


  for (int f = 0; f < nnotowned_bnd; f++) { CHECK(pr_id[f] >= 0); }
}
