/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"

#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"

#define SUPERMAP_TESTING 1
#include "SuperMapLumped.hh"

TEST(SUPERMAP_BLOCKMAP)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();
  int getSize = comm->getSize();

  if (getRank == 0)
    std::cout << "Test: SuperMapLumped with BlockMap with variable element size"
              << std::endl;

  // make a ghosted and local map 1
  std::vector<int> gids(3), size(3);
  for (int i = 0; i != 3; ++i) {
    gids[i] = 3 * getRank + i;
    size[i] = i + 1;
  }
  auto owned_map1 = Teuchos::rcp(
    new Epetra_BlockMap(3 * getSize, 3, &gids[0], &size[0], 0, *comm));
  CHECK_EQUAL(0, owned_map1->FirstPointInElement(0));
  CHECK_EQUAL(1, owned_map1->FirstPointInElement(1));
  CHECK_EQUAL(3, owned_map1->FirstPointInElement(2));

  if (getRank > 0) {
    gids.push_back(3 * getRank - 1);
    size.push_back(3);
  }
  if (getRank < getSize - 1) {
    gids.push_back(3 * (getRank + 1));
    size.push_back(1);
  }

  auto ghosted_map1 = Teuchos::rcp(
    new Epetra_BlockMap(-1, gids.size(), &gids[0], &size[0], 0, *comm));

  auto names = std::vector<std::string>{ "map1" };
  auto dofnums = std::vector<int>{ 1 };
  auto maps = std::vector<Teuchos::RCP<const Epetra_BlockMap>>{ owned_map1 };
  auto gmaps = std::vector<Teuchos::RCP<const Epetra_BlockMap>>{ ghosted_map1 };
  Operators::SuperMapLumped map(comm, names, dofnums, maps, gmaps);
}


TEST(SUPERMAP_BLOCKMAP_NULL_PROC)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();
  int getSize = comm->getSize();

  // only for parallel tests
  if (getSize > 1) {
    if (getRank == 0)
      std::cout
        << "Test: SuperMapLumped with BlockMap with variable element size"
        << std::endl;

    // make a ghosted and local map 1
    std::vector<int> gids(3), size(3);
    int n_cells;
    if (getRank != 1) {
      for (int i = 0; i != 3; ++i) {
        if (getRank == 0) {
          gids[i] = 3 * getRank + i;
          size[i] = i + 1;
        } else {
          gids[i] = 3 * (getRank - 1) + i;
          size[i] = i + 1;
        }
      }
      n_cells = 3;
    } else {
      n_cells = 0;
    }

    auto owned_map1 = Teuchos::rcp(new Epetra_BlockMap(
      3 * (getSize - 1), n_cells, &gids[0], &size[0], 0, *comm));
    if (getRank != 1) {
      CHECK_EQUAL(0, owned_map1->FirstPointInElement(0));
      CHECK_EQUAL(1, owned_map1->FirstPointInElement(1));
      CHECK_EQUAL(3, owned_map1->FirstPointInElement(2));
    } else {
      CHECK_EQUAL(-1, owned_map1->FirstPointInElement(0));
      CHECK_EQUAL(-1, owned_map1->FirstPointInElement(1));
      CHECK_EQUAL(-1, owned_map1->FirstPointInElement(2));
    }

    int n_cells_g = n_cells;
    if (getRank > 1) {
      gids.push_back(3 * (getRank - 1) - 1);
      size.push_back(3);
      n_cells_g++;
    }
    if ((getRank < getSize - 1)) {
      if (getRank > 1) {
        n_cells_g++;
        gids.push_back(3 * getRank);
        size.push_back(1);
      } else if (getRank == 0) {
        n_cells_g++;
        gids.push_back(3);
        size.push_back(1);
      }
    }

    auto ghosted_map1 =
      Teuchos::rcp(new Epetra_BlockMap(3 * (getSize - 1) + 2 * (getSize - 2),
                                       n_cells_g,
                                       &gids[0],
                                       &size[0],
                                       0,
                                       *comm));

    auto names = std::vector<std::string>{ "map1" };
    auto dofnums = std::vector<int>{ 1 };
    auto maps = std::vector<Teuchos::RCP<const Epetra_BlockMap>>{ owned_map1 };
    auto gmaps =
      std::vector<Teuchos::RCP<const Epetra_BlockMap>>{ ghosted_map1 };
    Operators::SuperMapLumped map(comm, names, dofnums, maps, gmaps);
  }
}


/* *****************************************************************
 * manually constructed test
 * *************************************************************** */
TEST(SUPERMAP_BLOCK_MANUAL)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();
  int getSize = comm->getSize();

  if (getRank == 0)
    std::cout << "Test: Manual test of SuperMapLumped with BLOCK maps"
              << std::endl;

  // make a ghosted and local map 1
  std::vector<int> gids(3), size(3);
  for (int i = 0; i != 3; ++i) {
    gids[i] = 3 * getRank + i;
    size[i] = i + 1;
  }

  auto owned_map1 = Teuchos::rcp(
    new Epetra_BlockMap(3 * getSize, 3, &gids[0], &size[0], 0, *comm));

  std::cout << *owned_map1 << "\n";
  std::cout << owned_map1->FirstPointInElement(0) << "\n";
  std::cout << owned_map1->FirstPointInElement(1) << "\n";
  std::cout << owned_map1->FirstPointInElement(2) << "\n";

  CHECK_EQUAL(0, owned_map1->FirstPointInElement(0));
  CHECK_EQUAL(1, owned_map1->FirstPointInElement(1));
  CHECK_EQUAL(3, owned_map1->FirstPointInElement(2));

  if (getRank > 0) {
    gids.push_back(3 * getRank - 1);
    size.push_back(3);
  }
  if (getRank < getSize - 1) {
    gids.push_back(3 * (getRank + 1));
    size.push_back(1);
  }

  auto ghosted_map1 = Teuchos::rcp(
    new Epetra_BlockMap(-1, gids.size(), &gids[0], &size[0], 0, *comm));

  // make a ghosted and local map 2
  std::vector<int> gids2(5), size2(5);
  for (int i = 0; i != 5; ++i) {
    gids2[i] = 5 * getRank + i;
    size2[i] = i % 2 + 1;
  }
  auto owned_map2 = Teuchos::rcp(
    new Epetra_BlockMap(5 * getSize, 5, &gids2[0], &size2[0], 0, *comm));

  if (getRank > 0) {
    gids2.push_back(5 * getRank - 1);
    size2.push_back(1);
  }
  if (getRank < getSize - 1) {
    gids2.push_back(5 * (getRank + 1));
    size2.push_back(1);
  }
  auto ghosted_map2 = Teuchos::rcp(
    new Epetra_BlockMap(-1, gids2.size(), &gids2[0], &size2[0], 0, *comm));

  // make the supermap
  std::vector<std::string> names;
  names.push_back("map1");
  names.push_back("map2");
  std::vector<int> dofnums(2, 1);

  std::vector<Teuchos::RCP<const Epetra_BlockMap>> maps;
  maps.push_back(owned_map1);
  maps.push_back(owned_map2);
  std::vector<Teuchos::RCP<const Epetra_BlockMap>> gmaps;
  gmaps.push_back(ghosted_map1);
  gmaps.push_back(ghosted_map2);

  Operators::SuperMapLumped map(comm, names, dofnums, maps, gmaps);

  std::cout << "======= Two Block Map =======" << std::endl;
  maps[0]->Print(std::cout);
  maps[1]->Print(std::cout);
  std::cout << "\n======= SuperMapLumped =======" << std::endl;
  map.getMap()->Print(std::cout);

  // check the offsets
  CHECK(map.Offset("map1") == 0);
  CHECK(map.Offset("map2") == 6);

  // check the indices
  {
    const std::vector<int>& inds_m1_d0 = map.Indices("map1", 0);
    // CHECK(inds_m1_d0.size() == 3);
    //    for (int i=0;i<3;i++) std::cout<<inds_m1_d0[i]<<" ";std::cout<<"\n";
    CHECK(inds_m1_d0[0] == 0);
    CHECK(inds_m1_d0[1] == 1);
    CHECK(inds_m1_d0[3] == 3);

    const std::vector<int>& inds_m2_d0 = map.Indices("map2", 0);
    // CHECK(inds_m2_d0.size() == 5);
    // for (int i=0;i<5;i++) std::cout<<inds_m2_d0[i]<<" ";std::cout<<"\n";
    CHECK(inds_m2_d0[0] == 6);
    CHECK(inds_m2_d0[1] == 7);
    CHECK(inds_m2_d0[3] == 9);
    CHECK(inds_m2_d0[4] == 10);
    CHECK(inds_m2_d0[6] == 12);
  }

  {
    const std::vector<int>& inds_m1_d0 = map.GhostIndices("map1", 0);
    CHECK(inds_m1_d0[0] == 0);
    CHECK(inds_m1_d0[1] == 1);
    CHECK(inds_m1_d0[3] == 3);

    const std::vector<int>& inds_m2_d0 = map.GhostIndices("map2", 0);
    // CHECK(inds_m2_d0.size() == 5);
    // for (int i=0;i<5;i++) std::cout<<inds_m2_d0[i]<<" ";std::cout<<"\n";
    CHECK(inds_m2_d0[0] == 6);
    CHECK(inds_m2_d0[1] == 7);
    CHECK(inds_m2_d0[3] == 9);
    CHECK(inds_m2_d0[4] == 10);
    CHECK(inds_m2_d0[6] == 12);
  }
}



