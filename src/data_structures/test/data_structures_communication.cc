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

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

using namespace Amanzi;

TEST(COMMUNICATION_PATTERN_DISTINCT_VECTORS)
{
  // having some trouble with communication -- lets make sure we understand
  // how to use tpetra
  auto comm = getDefaultComm();

  int size = comm->getSize();
  CHECK_EQUAL(2, size);
  int rank = comm->getRank();

  Map_ptr_type owned_map, ghost_map;
  if (rank == 0) {
    std::vector<int> gids_owned{ 0, 1 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 0, 1, 2 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  } else {
    std::vector<int> gids_owned{ 2, 3 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 2, 3, 1 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  }

  Import_type importer(owned_map, ghost_map);

  auto ghost = Teuchos::rcp(new MultiVector_type(ghost_map, 1));
  //  auto owned = ghost->offsetViewNonConst(owned_map, 0);
  auto owned = Teuchos::rcp(new MultiVector_type(owned_map, 1));

  owned->putScalar((double)(rank + 1));
  ghost->doImport(*owned, importer, Tpetra::INSERT);

  // manually get a host view and check
  ghost->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto ghost_v = ghost->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(2.0, ghost_v(2, 0), 1.e-6);
    } else {
      CHECK_CLOSE(2.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(1.0, ghost_v(2, 0), 1.e-6);
    }
  }
}


TEST(COMMUNICATION_PATTERN_OFFSET_VIEW)
{
  // having some trouble with communication -- lets make sure we understand
  // how to use tpetra
  auto comm = getDefaultComm();

  int size = comm->getSize();
  CHECK_EQUAL(2, size);
  int rank = comm->getRank();

  Map_ptr_type owned_map, ghost_map;
  if (rank == 0) {
    std::vector<int> gids_owned{ 0, 1 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 0, 1, 2 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  } else {
    std::vector<int> gids_owned{ 2, 3 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 2, 3, 1 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  }

  Import_type importer(owned_map, ghost_map);

  auto ghost = Teuchos::rcp(new MultiVector_type(ghost_map, 1));
  auto owned = ghost->offsetViewNonConst(owned_map, 0);

  owned->putScalar((double)(rank + 1));
  ghost->doImport(*owned, importer, Tpetra::INSERT);

  // manually get a host view and check
  ghost->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto ghost_v = ghost->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(2.0, ghost_v(2, 0), 1.e-6);
    } else {
      CHECK_CLOSE(2.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(1.0, ghost_v(2, 0), 1.e-6);
    }
  }
}


TEST(OFFSET_VIEW_1D_VIEW)
{
  // having some trouble with communication -- lets make sure we understand
  // how to use tpetra
  auto comm = getDefaultComm();

  int size = comm->getSize();
  CHECK_EQUAL(2, size);
  int rank = comm->getRank();
  int n_vecs = 5;

  Map_ptr_type owned_map, ghost_map;
  if (rank == 0) {
    std::vector<int> gids_owned{ 0, 1 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 0, 1, 2 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  } else {
    std::vector<int> gids_owned{ 2, 3 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 2, 3, 1 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  }

  auto ghost = Teuchos::rcp(new MultiVector_type(ghost_map, n_vecs));
  auto owned = ghost->offsetViewNonConst(owned_map, 0);

  CHECK_EQUAL(3, ghost->getLocalLength());
  CHECK_EQUAL(2, owned->getLocalLength());
  CHECK_EQUAL(n_vecs, ghost->getNumVectors());
  CHECK_EQUAL(n_vecs, owned->getNumVectors());

  {
    auto ghost_v = ghost->getLocalViewHost();
    CHECK_EQUAL(3, ghost_v.extent(0));
    CHECK_EQUAL(3 * n_vecs, ghost_v.span());
    CHECK_EQUAL(n_vecs, ghost_v.extent(1));

    auto owned_v = owned->getLocalViewHost();
    CHECK_EQUAL(2, owned_v.extent(0));

    // NOTE: this currently fails, and likely SHOULD fail according to Trilinos
    // CHECK_EQUAL(2*n_vecs, owned_v.span());
    CHECK_EQUAL(n_vecs, owned_v.extent(1));
  }

  CHECK(ghost->isConstantStride());
  CHECK_EQUAL(3 * n_vecs, ghost->get1dView().size());
  CHECK(owned->isConstantStride());
  // NOTE: this currently fails, and likely SHOULD NOT fail.  See Trilinos issue
  // #6058 CHECK_EQUAL(2*n_vecs, owned->get1dView().size());
}


TEST(COMMUNICATION_PATTERN_VANDELAY)
{
  // having some trouble with communication -- lets make sure we understand
  // how to use tpetra
  auto comm = getDefaultComm();

  int size = comm->getSize();
  CHECK_EQUAL(2, size);
  int rank = comm->getRank();

  Map_ptr_type owned_map, ghost_map, vandelay_map;
  if (rank == 0) {
    std::vector<int> gids_owned{ 0, 1 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 0, 1, 2 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));

    std::vector<int> gids_v{ 0 };
    vandelay_map = Teuchos::rcp(new Map_type(2, gids_v.data(), 1, 0, comm));

  } else {
    std::vector<int> gids_owned{ 2, 3 };
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{ 2, 3, 1 };
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));

    std::vector<int> gids_v{ 3 };
    vandelay_map = Teuchos::rcp(new Map_type(2, gids_v.data(), 1, 0, comm));
  }


  Import_type importer(owned_map, ghost_map);
  Import_type vandelay_importer(owned_map, vandelay_map);

  auto ghost = Teuchos::rcp(new MultiVector_type(ghost_map, 1));
  auto owned = ghost->offsetViewNonConst(owned_map, 0);
  auto vand = Teuchos::rcp(new MultiVector_type(vandelay_map, 1));

  owned->putScalar((double)(rank + 1));
  ghost->doImport(*owned, importer, Tpetra::INSERT);

  // manually get a host view and check
  ghost->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto ghost_v = ghost->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(2.0, ghost_v(2, 0), 1.e-6);
    } else {
      CHECK_CLOSE(2.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(1.0, ghost_v(2, 0), 1.e-6);
    }
  }

  // check the vandelay map
  vand->doImport(*owned, vandelay_importer, Tpetra::INSERT);
  vand->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto vand_v = vand->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, vand_v(0, 0), 1.e-6);
    } else {
      CHECK_CLOSE(2.0, vand_v(0, 0), 1.e-6);
    }
  }
}


TEST(COMMUNICATION_INTERMEDIATE)
{
  using namespace Amanzi;

  auto comm = getDefaultComm();
  int rank = comm->getRank();
  AmanziMesh::Preference pref;
  pref.clear();
  pref.push_back(AmanziMesh::Framework::MSTK);

  AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);

  auto mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 1, 1, 8);

  auto ghost = Teuchos::rcp(new MultiVector_type(mesh->face_map(true), 1));
  auto owned = ghost->offsetViewNonConst(mesh->face_map(false), 0);
  auto vand =
    Teuchos::rcp(new MultiVector_type(mesh->exterior_face_map(false), 1));

  Import_type importer(mesh->face_map(false), mesh->face_map(true));
  owned->putScalar((double)(rank + 1));
  ghost->doImport(*owned, importer, Tpetra::INSERT);

  // manually get a host view and check
  ghost->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto ghost_v = ghost->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(2.0,
                  ghost_v(mesh->num_entities(AmanziMesh::FACE,
                                             AmanziMesh::Parallel_type::OWNED),
                          0),
                  1.e-6);
    } else {
      CHECK_CLOSE(2.0, ghost_v(0, 0), 1.e-6);
      CHECK_CLOSE(1.0,
                  ghost_v(mesh->num_entities(AmanziMesh::FACE,
                                             AmanziMesh::Parallel_type::OWNED),
                          0),
                  1.e-6);
    }
  }

  // check the vandelay map
  Import_type vimporter(mesh->face_map(false), mesh->exterior_face_map(false));
  vand->doImport(*owned, vimporter, Tpetra::INSERT);
  // auto vimporter = mesh->exterior_face_importer();
  // vand->doImport(*owned, *vimporter, Tpetra::INSERT);
  vand->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto vand_v = vand->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, vand_v(0, 0), 1.e-6);
    } else {
      CHECK_CLOSE(2.0, vand_v(0, 0), 1.e-6);
    }
  }
}


TEST(IMPORTERS_BOUNDARY_FACES)
{
  // understand how to go both directions,
  //
  // bf = boundary_face_map.LID(face_map.GID(f));
  // and
  // f = face_map.LID(boundary_face_map.GID(bf));
  //
  // without calling doImport() (not always an option)
  std::cout << "IMPORTERS BOUNDARY FACES" << std::endl;  
  auto comm = getDefaultComm();

  int size = comm->getSize();
  CHECK_EQUAL(2, size);
  int rank = comm->getRank();

  Map_ptr_type boundary_face_map;

  auto face_map = Teuchos::rcp(new Map_type(8,4,0,comm));

  if (rank == 0) {
    std::vector<int> bf_gids{ 0, 3 };
    boundary_face_map = Teuchos::rcp(new Map_type(4, bf_gids.data(), 2, 0, comm));
  } else {
    std::vector<int> bf_gids{ 5, 6 };
    boundary_face_map = Teuchos::rcp(new Map_type(4, bf_gids.data(), 2, 0, comm));
  }

  Import_type importer(face_map, boundary_face_map);

  // go from bf to f
  IntVector_type vec(boundary_face_map);
  {
    auto vec_view = vec.getLocalViewDevice();
    auto import_view = importer.getPermuteFromLIDs_dv().view<DefaultHost>();
    auto import_to_view = importer.getPermuteToLIDs_dv().view<DefaultHost>();

    std::stringstream stream;
    stream << "Rank " << rank << ": len=" << import_view.extent(0) << "," << import_view.extent(1) << ": ";
    for (int i=0; i!=importer.getNumSameIDs(); ++i)
      stream << i << "-->" << i << ", ";
    for (int i=0; i!=import_view.extent(0); ++i)
      stream << import_view(i) << "-->" << import_to_view(i) << ", ";
    std::cout << stream.str() << std::endl;
    // assert(vec_view.extent(0) == import_view.extent(0));
    // Kokkos::parallel_for(vec_view.extent(0),
    //                      KOKKOS_LAMBDA(const int bf) {
    //                        vec_view(bf,0) = import_view(bf);
    //                      });
  }

  // {
  //   auto vec_view = vec.getLocalViewHost();
  //   for (int bf=0; bf!=vec_view.extent(0); ++bf) {
  //     std::cout << "bf: " << bf << " = " << vec_view(bf,0) << std::endl;
  //   }
  // }

}

