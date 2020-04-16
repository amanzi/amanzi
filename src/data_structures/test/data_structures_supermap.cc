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

#include "AmanziVector.hh"
#include "MeshFactory.hh"

#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector.hh"

#define SUPERMAP_TESTING 1
#include "SuperMapLumped.hh"

#include "SuperMap.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

Teuchos::RCP<Mesh>
getMesh(const Comm_ptr_type& comm)
{
  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list =
    plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  return meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
}


void
CHECK_OWNED_SUBSET_GHOST(cVectorView_type_<DefaultHost, int> owned,
                         int n_owned,
                         cVectorView_type_<DefaultHost, int> ghost,
                         int n_ghost)
{
  CHECK_EQUAL(n_owned, owned.size());
  CHECK_EQUAL(n_ghost, ghost.size());
  for (int i = 0; i != n_owned; ++i) { CHECK_EQUAL(owned[i], ghost[i]); }
}

void
CHECK_UNIQUE(std::vector<cVectorView_type_<DefaultHost, int>> index_lists,
             int total)
{
  std::set<int> all;
  int expected_count = 0;
  for (auto index_list : index_lists) {
    expected_count += index_list.size();
    for (std::size_t i = 0; i != index_list.extent(0); ++i)
      all.emplace(index_list(i));
  }
  CHECK_EQUAL(expected_count, all.size());
  CHECK_EQUAL(total, expected_count);
}

void
CHECK_SAME(cVectorView_type_<DefaultHost, int> one,
           cVectorView_type_<DefaultHost, int> two)
{
  CHECK_EQUAL(one.extent(0), two.extent(0));
  for (std::size_t i = 0; i != one.extent(0); ++i) {
    CHECK_EQUAL(one[i], two[i]);
  }
}


/* *****************************************************************
 * manually constructed test
 * **************************************************************** */
TEST(SUPERMAP_MANUAL)
{
  //  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();
  int getSize = comm->getSize();

  if (getRank == 0)
    std::cout << "Test: Manual test of SuperMapLumped" << std::endl;

  // make a ghosted and local map 1
  std::vector<int> gids(3);
  for (int i = 0; i != 3; ++i) { gids[i] = 3 * getRank + i; }
  auto owned_map1 =
    Teuchos::rcp(new Map_type(3 * getSize, gids.data(), 3, 0, comm));

  if (getRank > 0) gids.push_back(3 * getRank - 1);
  if (getRank < getSize - 1) gids.push_back(3 * (getRank + 1));
  auto ghosted_map1 =
    Teuchos::rcp(new Map_type(-1, gids.data(), gids.size(), 0, comm));

  // make a ghosted and local map 2
  std::vector<int> gids2(5);
  for (int i = 0; i != 5; ++i) { gids2[i] = 5 * getRank + i; }
  auto owned_map2 =
    Teuchos::rcp(new Map_type(5 * getSize, gids2.data(), 5, 0, comm));

  if (getRank > 0) gids2.push_back(5 * getRank - 1);
  if (getRank < getSize - 1) gids2.push_back(5 * (getRank + 1));
  auto ghosted_map2 =
    Teuchos::rcp(new Map_type(-1, gids2.data(), gids2.size(), 0, comm));

  // make the supermap
  std::vector<std::string> names = { "map1", "map2" };
  std::map<std::string, std::size_t> dofnums;
  dofnums["map1"] = 2;
  dofnums["map2"] = 2;
  std::map<std::string, Teuchos::RCP<const BlockMap_type>> maps;
  maps["map1"] = owned_map1;
  maps["map2"] = owned_map2;

  std::map<std::string, Teuchos::RCP<const BlockMap_type>> gmaps;
  gmaps["map1"] = ghosted_map1;
  gmaps["map2"] = ghosted_map2;

  Operators::SuperMapLumped map(
    Teuchos::rcp(new BlockSpace(comm, names, maps, gmaps, dofnums)));

  // check the offsets
  CHECK(map.Offset("map1") == 0);
  CHECK(map.Offset("map2") == 2 * 3);

  // check the ghosted offsets
  CHECK(map.GhostedOffset("map1") == (2 * 3 + 2 * 5));
  if (getSize > 1) {
    if (getRank == 0 || getRank == getSize - 1) {
      CHECK(map.GhostedOffset("map2") == (2 * 3 + 2 * 5 + 2));
    } else {
      CHECK(map.GhostedOffset("map2") == (2 * 3 + 2 * 5 + 4));
    }
  } else {
    CHECK(map.GhostedOffset("map2") == (2 * 3 + 2 * 5));
  }

  // check num owned
  CHECK(map.NumOwnedElements("map1") == 3);
  CHECK(map.NumOwnedElements("map2") == 5);

  // check num owned/used
  if (getSize > 1) {
    if (getRank == 0 || getRank == getSize - 1) {
      CHECK(map.NumUsedElements("map1") == 4);
      CHECK(map.NumUsedElements("map2") == 6);
    } else {
      CHECK(map.NumUsedElements("map1") == 5);
      CHECK(map.NumUsedElements("map2") == 7);
    }
  } else {
    CHECK(map.NumUsedElements("map1") == 3);
    CHECK(map.NumUsedElements("map2") == 5);
  }

  // check num dofs
  CHECK(map.getNumVectors("map1") == 2);
  CHECK(map.getNumVectors("map2") == 2);

  // check that the maps properly export
  Vector_type_<int> ghosted(map.getGhostedMap());
  Vector_type_<int> owned(map.getMap());
  Import_type importer(map.getMap(), map.getGhostedMap());

  {
    auto owned_v = owned.getLocalViewHost();
    for (int i = 0; i != map.getMap()->getNodeNumElements(); ++i) {
      owned_v(i, 0) = map.getMap()->getGlobalElement(i);
    }
  }

  ghosted.putScalar(0);
  ghosted.doImport(owned, importer, Tpetra::INSERT);
  ghosted.sync_host();
  {
    auto ghosted_v = ghosted.getLocalViewHost();
    for (int i = 0; i != map.getGhostedMap()->getNodeNumElements(); ++i) {
      CHECK_EQUAL(ghosted_v(i, 0), map.getGhostedMap()->getGlobalElement(i));
    }
  }

  // check the indices
  {
    auto inds_m1_d0 = map.Indices<DefaultHost>("map1", 0);
    CHECK_EQUAL(3, inds_m1_d0.size());
    CHECK_EQUAL(0, inds_m1_d0[0]);
    CHECK_EQUAL(2, inds_m1_d0[1]);
    CHECK_EQUAL(4, inds_m1_d0[2]);

    auto inds_m1_d1 = map.Indices<DefaultHost>("map1", 1);
    CHECK_EQUAL(3, inds_m1_d1.size());
    CHECK_EQUAL(1, inds_m1_d1[0]);
    CHECK_EQUAL(3, inds_m1_d1[1]);
    CHECK_EQUAL(5, inds_m1_d1[2]);

    auto inds_m2_d0 = map.Indices<DefaultHost>("map2", 0);
    CHECK_EQUAL(5, inds_m2_d0.size());
    CHECK_EQUAL(6, inds_m2_d0[0]);
    CHECK_EQUAL(8, inds_m2_d0[1]);
    CHECK_EQUAL(10, inds_m2_d0[2]);
    CHECK_EQUAL(12, inds_m2_d0[3]);
    CHECK_EQUAL(14, inds_m2_d0[4]);

    auto inds_m2_d1 = map.Indices<DefaultHost>("map2", 1);
    CHECK(inds_m2_d1.size() == 5);
    CHECK_EQUAL(7, inds_m2_d1[0]);
    CHECK_EQUAL(9, inds_m2_d1[1]);
    CHECK_EQUAL(11, inds_m2_d1[2]);
    CHECK_EQUAL(13, inds_m2_d1[3]);
    CHECK_EQUAL(15, inds_m2_d1[4]);
  }

  {
    auto inds_m1_d0 = map.GhostIndices<DefaultHost>("map1", 0);
    CHECK_EQUAL(0, inds_m1_d0[0]);
    CHECK_EQUAL(2, inds_m1_d0[1]);
    CHECK_EQUAL(4, inds_m1_d0[2]);

    auto inds_m1_d1 = map.GhostIndices<DefaultHost>("map1", 1);
    CHECK_EQUAL(1, inds_m1_d1[0]);
    CHECK_EQUAL(3, inds_m1_d1[1]);
    CHECK_EQUAL(5, inds_m1_d1[2]);

    auto inds_m2_d0 = map.GhostIndices<DefaultHost>("map2", 0);
    CHECK_EQUAL(6, inds_m2_d0[0]);
    CHECK_EQUAL(8, inds_m2_d0[1]);
    CHECK_EQUAL(10, inds_m2_d0[2]);
    CHECK_EQUAL(12, inds_m2_d0[3]);
    CHECK_EQUAL(14, inds_m2_d0[4]);

    auto inds_m2_d1 = map.GhostIndices<DefaultHost>("map2", 1);
    CHECK_EQUAL(7, inds_m2_d1[0]);
    CHECK_EQUAL(9, inds_m2_d1[1]);
    CHECK_EQUAL(11, inds_m2_d1[2]);
    CHECK_EQUAL(13, inds_m2_d1[3]);
    CHECK_EQUAL(15, inds_m2_d1[4]);

    if (getSize > 1) {
      if (getRank == 0 || getRank == getSize - 1) {
        CHECK_EQUAL(4, inds_m1_d0.size());
        CHECK_EQUAL(4, inds_m1_d1.size());
        CHECK_EQUAL(6, inds_m2_d0.size());
        CHECK_EQUAL(6, inds_m2_d1.size());

        CHECK_EQUAL(16, inds_m1_d0[3]);
        CHECK_EQUAL(17, inds_m1_d1[3]);
        CHECK_EQUAL(18, inds_m2_d0[5]);
        CHECK_EQUAL(19, inds_m2_d1[5]);
      } else {
        CHECK_EQUAL(5, inds_m1_d0.size());
        CHECK_EQUAL(5, inds_m1_d1.size());
        CHECK_EQUAL(7, inds_m2_d0.size());
        CHECK_EQUAL(7, inds_m2_d1.size());

        CHECK_EQUAL(16, inds_m1_d0[3]);
        CHECK_EQUAL(18, inds_m1_d0[4]);
        CHECK_EQUAL(17, inds_m1_d1[3]);
        CHECK_EQUAL(19, inds_m1_d1[4]);
        CHECK_EQUAL(20, inds_m2_d0[5]);
        CHECK_EQUAL(22, inds_m2_d0[6]);
        CHECK_EQUAL(21, inds_m2_d1[5]);
        CHECK_EQUAL(23, inds_m2_d1[6]);
      }
    } else {
      CHECK_EQUAL(3, inds_m1_d0.size());
      CHECK_EQUAL(3, inds_m1_d1.size());
      CHECK_EQUAL(5, inds_m2_d0.size());
      CHECK_EQUAL(5, inds_m2_d1.size());
    }
  }
}


TEST(SUPERMAP_FROM_SINGLE_COMPOSITEVECTOR)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();
  if (getRank == 0)
    std::cout << "Test: SuperMapLumped from 1 CompositeVector with multiple "
                 "components and multiple dofs"
              << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cv;
  std::vector<std::string> names = { "cell", "face" };
  std::vector<int> dofs{ 2, 2 };
  std::vector<Entity_kind> locs = { Entity_kind::CELL, Entity_kind::FACE };
  cv.SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map =
    createSuperMap(cv.CreateSpace().ptr());

  int ncells_owned = mesh->num_entities(CELL, Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(FACE, Parallel_type::OWNED);
  int ncells_used = mesh->num_entities(CELL, Parallel_type::ALL);
  int nfaces_used = mesh->num_entities(FACE, Parallel_type::ALL);

  // check basic sizes
  CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned,
              map->getMap()->getNodeNumElements());
  CHECK_EQUAL(2 * ncells_used + 2 * nfaces_used,
              map->getGhostedMap()->getNodeNumElements());

  // check CompMaps
  CHECK(mesh->cell_map(false)->isSameAs(*map->ComponentMap(0, "cell")));
  CHECK(mesh->cell_map(true)->isSameAs(*map->ComponentGhostedMap(0, "cell")));
  CHECK(mesh->face_map(false)->isSameAs(*map->ComponentMap(0, "face")));
  CHECK(mesh->face_map(true)->isSameAs(*map->ComponentGhostedMap(0, "face")));

  // check ordering is as expected
  const auto& inds_c0 = map->GhostIndices<DefaultHost>(0, "cell", 0);
  const auto& inds_c1 = map->GhostIndices<DefaultHost>(0, "cell", 1);
  const auto& inds_f0 = map->GhostIndices<DefaultHost>(0, "face", 0);
  const auto& inds_f1 = map->GhostIndices<DefaultHost>(0, "face", 1);

  const auto& inds_c0_owned = map->Indices<DefaultHost>(0, "cell", 0);
  const auto& inds_c1_owned = map->Indices<DefaultHost>(0, "cell", 1);
  const auto& inds_f0_owned = map->Indices<DefaultHost>(0, "face", 0);
  const auto& inds_f1_owned = map->Indices<DefaultHost>(0, "face", 1);

  // check owned list is a subset of the ghost list
  CHECK_OWNED_SUBSET_GHOST(inds_c0_owned, ncells_owned, inds_c0, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_c1_owned, ncells_owned, inds_c1, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_f0_owned, nfaces_owned, inds_f0, nfaces_used);
  CHECK_OWNED_SUBSET_GHOST(inds_f1_owned, nfaces_owned, inds_f1, nfaces_used);

  // check lists are non-overlapping and of sufficient size to cover all values
  CHECK_UNIQUE(
    std::vector<cVectorView_type_<DefaultHost, int>>{
      inds_c0_owned, inds_c1_owned, inds_f0_owned, inds_f1_owned },
    2 * ncells_owned + 2 * nfaces_owned);
  CHECK_UNIQUE(
    std::vector<cVectorView_type_<DefaultHost, int>>{
      inds_c0, inds_c1, inds_f0, inds_f1 },
    2 * ncells_used + 2 * nfaces_used);

  // check owned maps interleave
  CHECK_EQUAL(0, inds_c0[0]);
  CHECK_EQUAL(1, inds_c1[0]);
  CHECK_EQUAL(2, inds_c0[1]);
  CHECK_EQUAL(3, inds_c1[1]);

  CHECK_EQUAL(2 * ncells_owned, inds_f0[0]);
  CHECK_EQUAL(2 * ncells_owned + 1, inds_f1[0]);
  CHECK_EQUAL(2 * ncells_owned + 2, inds_f0[1]);
  CHECK_EQUAL(2 * ncells_owned + 3, inds_f1[1]);


  // check ghosts pick up at the end of owned and interleave
  if (comm->getSize() > 1) {
    CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned, inds_c0[ncells_owned]);
    CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned + 1, inds_c1[ncells_owned]);
    CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned + 2,
                inds_c0[ncells_owned + 1]);
    CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned + 3,
                inds_c1[ncells_owned + 1]);

    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned, inds_f0[nfaces_owned]);
    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned + 1, inds_f1[nfaces_owned]);
    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned + 2,
                inds_f0[nfaces_owned + 1]);
    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned + 3,
                inds_f1[nfaces_owned + 1]);
  }
}


TEST(SUPERMAP_FROM_SINGLE_COMPOSITEVECTOR_REPEATED_MAPS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();
  if (getRank == 0)
    std::cout << "Test: SuperMapLumped from 1 CompositeVector with multiple "
                 "components which include repeated maps, should still "
                 "interleave as if this were one map and two dofs"
              << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cv;
  std::vector<std::string> names;
  names.push_back("cellA");
  names.push_back("cellB");
  std::vector<int> dofs{ 1, 1 };
  std::vector<Entity_kind> locs;
  locs.push_back(CELL);
  locs.push_back(CELL);
  cv.SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map =
    createSuperMap(cv.CreateSpace().ptr());

  int ncells_owned = mesh->num_entities(CELL, Parallel_type::OWNED);
  int ncells_used = mesh->num_entities(CELL, Parallel_type::ALL);

  // check basic sizes
  CHECK_EQUAL(2 * ncells_owned, map->getMap()->getNodeNumElements());
  CHECK_EQUAL(2 * ncells_used, map->getGhostedMap()->getNodeNumElements());

  // check CompMaps
  CHECK(mesh->cell_map(false)->isSameAs(*map->ComponentMap(0, "cellA")));
  CHECK(mesh->cell_map(true)->isSameAs(*map->ComponentGhostedMap(0, "cellA")));
  CHECK(mesh->cell_map(false)->isSameAs(*map->ComponentMap(0, "cellB")));
  CHECK(mesh->cell_map(true)->isSameAs(*map->ComponentGhostedMap(0, "cellB")));

  // check ordering is as expected
  const auto& inds_c0 = map->GhostIndices<DefaultHost>(0, "cellA", 0);
  const auto& inds_c1 = map->GhostIndices<DefaultHost>(0, "cellB", 0);

  const auto& inds_c0_owned = map->Indices<DefaultHost>(0, "cellA", 0);
  const auto& inds_c1_owned = map->Indices<DefaultHost>(0, "cellB", 0);

  // check owned list is a subset of the ghost list
  CHECK_OWNED_SUBSET_GHOST(inds_c0_owned, ncells_owned, inds_c0, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_c1_owned, ncells_owned, inds_c1, ncells_used);

  // check lists are non-overlapping and of sufficient size to cover all values
  CHECK_UNIQUE(
    std::vector<cVectorView_type_<DefaultHost, int>>{ inds_c0_owned,
                                                            inds_c1_owned },
    2 * ncells_owned);
  CHECK_UNIQUE(
    std::vector<cVectorView_type_<DefaultHost, int>>{ inds_c0, inds_c1 },
    2 * ncells_used);

  // check owned maps interleave
  CHECK_EQUAL(0, inds_c0[0]);
  CHECK_EQUAL(1, inds_c1[0]);
  CHECK_EQUAL(2, inds_c0[1]);
  CHECK_EQUAL(3, inds_c1[1]);

  // check ghosts pick up at the end of owned and interleave
  if (comm->getSize() > 1) {
    CHECK_EQUAL(2 * ncells_owned, inds_c0[ncells_owned]);
    CHECK_EQUAL(2 * ncells_owned + 1, inds_c1[ncells_owned]);
    CHECK_EQUAL(2 * ncells_owned + 2, inds_c0[ncells_owned + 1]);
    CHECK_EQUAL(2 * ncells_owned + 3, inds_c1[ncells_owned + 1]);
  }

  // this map should be IDENTICAL to that of single component, 2 dofs
  CompositeVectorSpace cv2;
  std::vector<std::string> names2{ "cell" };
  std::vector<int> dofs2{ 2 };
  std::vector<Entity_kind> locs2{ CELL };
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents(names2, locs2, dofs2);

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map2 =
    createSuperMap(cv2.CreateSpace().ptr());

  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));
}


TEST(SUPERMAP_FROM_TWO_IDENTICAL_COMPOSITEVECTORS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();
  if (getRank == 0)
    std::cout << "Test: SuperMapLumped from 2 CompositeVectors with same map, "
                 "single dof is same as 1 CompositeVector with two dofs"
              << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 1 });
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 1 });

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(
    new SuperMap(comm, { cvA.CreateSpace().ptr(), cvB.CreateSpace().ptr() }));

  // now create another with one CV, 2 dofs
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 2 });
  Teuchos::RCP<Operators::SuperMap> map2 =
    createSuperMap(cv2.CreateSpace().ptr());

  // same map!
  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));

  // same indices!
  CHECK_SAME(map->Indices<DefaultHost>(0, "cell", 0),
             map2->Indices<DefaultHost>(0, "cell", 0));
  CHECK_SAME(map->GhostIndices<DefaultHost>(0, "cell", 0),
             map2->GhostIndices<DefaultHost>(0, "cell", 0));
  CHECK_SAME(map->Indices<DefaultHost>(1, "cell", 0),
             map2->Indices<DefaultHost>(0, "cell", 1));
  CHECK_SAME(map->GhostIndices<DefaultHost>(1, "cell", 0),
             map2->GhostIndices<DefaultHost>(0, "cell", 1));
}


TEST(SUPERMAP_FROM_CELL_PLUS_FACE_IS_CELLFACE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();
  if (getRank == 0)
    std::cout
      << "Test: SuperMapLumped from 2 CompositeVectors with different maps, "
         "single dof is same as 1 CompositeVector with two components"
      << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell" },
    std::vector<Entity_kind>{ CELL },
    std::vector<int>{ 1 });
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "face" },
    std::vector<Entity_kind>{ FACE },
    std::vector<int>{ 1 });
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell", "face" },
    std::vector<Entity_kind>{ CELL, FACE },
    std::vector<int>{ 1, 1 });

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(
    new SuperMap(comm, { cvA.CreateSpace().ptr(), cvB.CreateSpace().ptr() }));
  Teuchos::RCP<Operators::SuperMap> map2 =
    Teuchos::rcp(new SuperMap(comm, { cv2.CreateSpace().ptr() }));

  // same map!
  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));

  // same indices!
  CHECK_SAME(map->Indices<DefaultHost>(0, "cell", 0),
             map2->Indices<DefaultHost>(0, "cell", 0));
  CHECK_SAME(map->GhostIndices<DefaultHost>(0, "cell", 0),
             map2->GhostIndices<DefaultHost>(0, "cell", 0));
  CHECK_SAME(map->Indices<DefaultHost>(1, "face", 0),
             map2->Indices<DefaultHost>(0, "face", 0));
  CHECK_SAME(map->GhostIndices<DefaultHost>(1, "face", 0),
             map2->GhostIndices<DefaultHost>(0, "face", 0));
}


TEST(SUPERMAP_FROM_SAME_NAME_DIFFERENT_MAP)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();
  if (getRank == 0)
    std::cout
      << "Test: SuperMapLumped from 2 CompositeVectors with different maps but "
         "the same name, is same as 1 CompositeVector with two components"
      << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell" },
    std::vector<Entity_kind>{ CELL },
    std::vector<int>{ 1 });
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell" },
    std::vector<Entity_kind>{ FACE },
    std::vector<int>{ 1 });
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell", "face" },
    std::vector<Entity_kind>{ CELL, FACE },
    std::vector<int>{ 1, 1 });

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(
    new SuperMap(comm, { cvA.CreateSpace().ptr(), cvB.CreateSpace().ptr() }));
  Teuchos::RCP<Operators::SuperMap> map2 =
    Teuchos::rcp(new SuperMap(comm, { cv2.CreateSpace().ptr() }));

  // same map!
  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));

  // same indices!
  CHECK_SAME(map->Indices<DefaultHost>(0, "cell", 0),
             map2->Indices<DefaultHost>(0, "cell", 0));
  CHECK_SAME(map->GhostIndices<DefaultHost>(0, "cell", 0),
             map2->GhostIndices<DefaultHost>(0, "cell", 0));
  CHECK_SAME(map->Indices<DefaultHost>(1, "cell", 0),
             map2->Indices<DefaultHost>(0, "face", 0));
  CHECK_SAME(map->GhostIndices<DefaultHost>(1, "cell", 0),
             map2->GhostIndices<DefaultHost>(0, "face", 0));
}


//
// these are block tests, and cannot pass as there is no Tpetra::BlockMap

// TEST(SUPERMAP_FROM_SAME_NAME_SAME_MAP_DIFFERENT_ELEMENTSIZE) {
//   //
//   // This test should be good for pressure in matrix + fracture?
//   //
//   using namespace Amanzi;
//   using namespace Amanzi::AmanziMesh;
//   using namespace Amanzi::AmanziGeometry;
//   using namespace Amanzi::Operators;

//   auto comm = getDefaultComm();
//   int getRank = comm->getRank();
//   if (getRank == 0) std::cout << "Test: SuperMapLumped from 2
//   CompositeVectors with same elements but different element sizes in the same
//   compname" << std::endl;

//   auto mesh = getMesh(comm);
//   int ncells = mesh->num_entities(CELL, Parallel_type::OWNED);

//   // create a CVSpace
//   const auto& cell_map = mesh->cell_map(false);
//   const int * gids = nullptr;
//   const long long * llgids = nullptr;
//   cell_map.MyGlobalElements(gids, llgids);
//   Epetra_IntVector element_size(cell_map);
//   element_size.PutValue(1);
//   element_size[1] = 2;

//   Teuchos::RCP<const Epetra_BlockMap> block_cell_map = Teuchos::rcp(new
//   Epetra_BlockMap(cell_map.NumGlobalElements(),
//   cell_map.getNodeNumElements(), gids, &element_size[0], 0, comm));

//   const auto& cell_mapg = mesh->cell_map(true);
//   Epetra_Import importer(cell_map, cell_mapg);
//   Epetra_IntVector element_sizeg(cell_mapg);
//   element_sizeg.Import(element_size, importer, Insert);

//   cell_mapg.MyGlobalElements(gids, llgids);
//   Teuchos::RCP<const Epetra_BlockMap> block_cell_map_g = Teuchos::rcp(new
//   Epetra_BlockMap(cell_mapg.NumGlobalElements(),
//   cell_mapg.getNodeNumElements(), gids, &element_sizeg[0], 0, comm));


//   CompositeVectorSpace cvA;
//   cvA.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL}, {1});
//   CompositeVectorSpace cvB;
//   cvB.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL},
//   {{"cell",block_cell_map}}, {{"cell", block_cell_map_g}}, {1});

//   // create a SuperMapLumped from this space
//   Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(new SuperMap(comm,
//   {cvA.CreateSpace().ptr(), cvB.CreateSpace().ptr()}));

//   // not the same map!
//   CHECK(map->Indices<DefaultHost>(0, "cell", 0).size() == ncells);
//   CHECK(map->Indices<DefaultHost>(1, "cell", 0).size() ==
//   block_cell_map->NumMyPoints()); CHECK_EQUAL(1 + ncells,
//   block_cell_map->NumMyPoints());

//   const auto& inds1 = map->Indices<DefaultHost>(0, "cell", 0);
//   const auto& inds2 = map->Indices<DefaultHost>(1, "cell", 0);
//   CHECK_EQUAL(0, inds1[0]);
//   CHECK_EQUAL(1, inds1[1]);

//   CHECK_EQUAL(ncells, inds2[0]);
//   CHECK_EQUAL(ncells+1, inds2[1]);

//   CHECK_UNIQUE(std::vector<cVectorView_type_<DefaultHost,int> >{ inds1,
//   inds2 }, 2*ncells + 1);

// }


TEST(SUPERMAP_FROM_TREEVECTOR)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "Test: FD like matrix, null off-proc assembly" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list =
    plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh =
  //  meshfactory.create("test/median32x33.exo");

  // create a CVSpace
  Teuchos::RCP<CompositeVectorSpace> cv =
    Teuchos::rcp(new CompositeVectorSpace());
  std::vector<std::string> names;
  names.push_back("cell");
  names.push_back("face");
  std::vector<int> dofs(2, 2);
  std::vector<Entity_kind> locs;
  locs.push_back(CELL);
  locs.push_back(FACE);
  cv->SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a TV
  TreeVectorSpace tv;
  Teuchos::RCP<TreeVectorSpace> tv_ss1 = Teuchos::rcp(new TreeVectorSpace());
  tv_ss1->SetData(cv->CreateSpace());
  Teuchos::RCP<TreeVectorSpace> tv_ss2 = Teuchos::rcp(new TreeVectorSpace());
  tv_ss2->SetData(cv->CreateSpace());
  tv.PushBack(tv_ss1);
  tv.PushBack(tv_ss2);

  // create a SuperMapLumped from a singleton space
  Teuchos::RCP<Operators::SuperMap> map_singleton = createSuperMap(*tv_ss1);
  // create a SuperMapLumped from a tree space
  Teuchos::RCP<Operators::SuperMap> map = createSuperMap(tv);
}


//
// Tests copying to and from other vectors
//
struct Maps {
  Maps()
  {
    comm = Amanzi::getDefaultComm();

    // create a mesh
    mesh = getMesh(comm);

    // create a vector
    cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1);

    Teuchos::RCP<TreeVectorSpace> tvs0 = Teuchos::rcp(new TreeVectorSpace());
    tvs0->SetData(cvs->CreateSpace());

    tvs = Teuchos::rcp(new TreeVectorSpace());
    tvs->PushBack(tvs0);
    tvs->PushBack(tvs0);

    // create a supermap, vec
    map = createSuperMap(*tvs);
  }

  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<CompositeVectorSpace> cvs;
  Teuchos::RCP<TreeVectorSpace> tvs;
  Teuchos::RCP<SuperMap> map;
};


TEST(SUPERMAP_COPY_INVERTIBLE)
{
  Maps maps;
  Teuchos::RCP<TreeVector> tv = Teuchos::rcp(new TreeVector(maps.tvs));

  // initialize randomly
  tv->SubVector(0)->Data()->random();
  tv->SubVector(1)->Data()->random();

  Vector_type vec(maps.map->getMap());

  // copy forward, backward
  Teuchos::RCP<TreeVector> tv2 = Teuchos::rcp(new TreeVector(maps.tvs));
  int ierr = copyToSuperVector(*maps.map, *tv, vec);
  CHECK(!ierr);

  ierr = copyFromSuperVector(*maps.map, vec, *tv2);
  CHECK(!ierr);

  // check the same
  tv2->update(-1., *tv, 1.);
  double norm;
  norm = tv2->norm2();
  CHECK_CLOSE(0., norm, 1.e-16);
}


TEST(SUPERMAP_COPY_INTS)
{
  Maps maps;
  Teuchos::RCP<TreeVector> tv = Teuchos::rcp(new TreeVector(maps.tvs));
  tv->SubVector(0)->Data()->GetComponent("cell", false)->putScalar(3.);
  tv->SubVector(1)->Data()->GetComponent("cell", false)->putScalar(4.);
  tv->SubVector(0)->Data()->GetComponent("face", false)->putScalar(5.);
  tv->SubVector(1)->Data()->GetComponent("face", false)->putScalar(6.);

  Vector_type vec(maps.map->getMap());

  // copy forward
  Teuchos::RCP<TreeVector> tv2 = Teuchos::rcp(new TreeVector(maps.tvs));
  int ierr = copyToSuperVector(*maps.map, *tv, vec);
  CHECK(!ierr);

  // check values
  int ncells =
    maps.mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces =
    maps.mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // check sizes
  CHECK_EQUAL(2 * ncells + 2 * nfaces, vec.getLocalLength());

  {
    auto vec_v = vec.getLocalViewHost();
    for (int i = 0; i != ncells; ++i) {
      CHECK_EQUAL(3., vec_v(i * 2,0));
      CHECK_EQUAL(4., vec_v(i * 2 + 1,0));
    }

    for (int i = 0; i != nfaces; ++i) {
      CHECK_EQUAL(5., vec_v(2 * ncells + i * 2,0));
      CHECK_EQUAL(6., vec_v(2 * ncells + i * 2 + 1,0));
    }
  }
}
