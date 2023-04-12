/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

*/

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
#include "Iterators.hh"

#define SUPERMAP_TESTING 1
#include "SuperMapLumped.hh"
#include "SuperMap.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;


Teuchos::RCP<Mesh>
getMesh(const Comm_ptr_type& comm)
{
  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  return meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
}


template<class View_type, typename size_type>
void
CHECK_OWNED_SUBSET_GHOST(const View_type& owned,
                         size_type n_owned,
                         const View_type& ghost,
                         size_type n_ghost)
{
  CHECK(n_owned == owned.extent(0));
  CHECK(n_ghost == ghost.extent(0));
  for (int i = 0; i != n_owned; ++i) { CHECK_EQUAL(owned(i), ghost(i)); }
}

template<class View_type, typename size_type>
void
CHECK_UNIQUE(const std::vector<View_type>& index_lists, size_type total)
{
  std::set<Entity_GID> all;
  int expected_count = 0;
  for (auto index_list : index_lists) {
    expected_count += index_list->extent(0);
    all.insert(begin(*index_list), end(*index_list));
  }
  CHECK_EQUAL(expected_count, all.size());
  CHECK_EQUAL(total, expected_count);
}

template<class View_type>
bool
CHECK_EQUAL_VIEW(const View_type& v1, const View_type& v2) {
  if (v1.size() != v2.size()) return false;
  auto v1i = begin(v1);
  auto v2i = begin(v2);
  for (int i=0; i!=v1.size(); ++i) {
    if (*v1i != *v2i) return false;
    v1i++;
    v2i++;
  }
  return true;
}


/* *****************************************************************
 * manually constructed test
 * **************************************************************** */
void
SuperMap_Manual(bool continuous)
{
  auto comm = getDefaultComm();
  int MyPID = comm->getRank();
  int NumProc = comm->getSize();

  if (MyPID == 0) std::cout << "Test: Manual test of SuperMapLumped" << std::endl;

  // make a ghosted and local map 1
  std::vector<int> gids(3);
  if (continuous) {
    for (int i = 0; i != 3; ++i) { gids[i] = 3 * MyPID + i; }
  } else {
    for (int i = 0; i != 3; ++i) { gids[i] = 3 * (NumProc - MyPID) + i; }
  }
  Map_ptr_type owned_map1 =
    Teuchos::rcp(new Map_type(3 * NumProc, gids, 0, comm));

  if (MyPID > 0) gids.push_back(3 * MyPID - 1);
  if (MyPID < NumProc - 1) gids.push_back(3 * (MyPID + 1));
  Map_ptr_type ghosted_map1 =
    Teuchos::rcp(new Map_type(-1, gids, 0, comm));

  // make a ghosted and local map 2
  std::vector<int> gids2(5);
  for (int i = 0; i != 5; ++i) { gids2[i] = 5 * MyPID + i; }
  Map_ptr_type owned_map2 =
    Teuchos::rcp(new Map_type(5 * NumProc, gids2, 0, comm));

  if (MyPID > 0) gids2.push_back(5 * MyPID - 1);
  if (MyPID < NumProc - 1) gids2.push_back(5 * (MyPID + 1));
  Map_ptr_type ghosted_map2 =
    Teuchos::rcp(new Map_type(-1, gids2, 0, comm));

  // make the supermap
  std::vector<std::string> names;
  names.push_back("map1");
  names.push_back("map2");
  std::map<std::string,std::size_t> dofnums;
  dofnums["map1"] = 2;
  dofnums["map2"] = 2;

  std::map<std::string, BlockMap_ptr_type> maps;
  maps["map1"] = owned_map1;
  maps["map2"] = owned_map2;
  std::map<std::string, BlockMap_ptr_type> gmaps;
  gmaps["map1"] = ghosted_map1;
  gmaps["map2"] = ghosted_map2;

  auto space = Teuchos::rcp(new BlockSpace(comm, names, maps, gmaps, dofnums));
  SuperMapLumped map(space);

  // check the offsets
  CHECK(map.getOffset("map1") == 0);
  CHECK(map.getOffset("map2") == 2 * 3);

  // check the ghosted offsets
  CHECK(map.getGhostedOffset("map1") == (2 * 3 + 2 * 5));
  if (NumProc > 1) {
    if (MyPID == 0 || MyPID == NumProc - 1) {
      CHECK(map.getGhostedOffset("map2") == (2 * 3 + 2 * 5 + 2));
    } else {
      CHECK(map.getGhostedOffset("map2") == (2 * 3 + 2 * 5 + 4));
    }
  } else {
    CHECK(map.getGhostedOffset("map2") == (2 * 3 + 2 * 5));
  }

  // check num owned
  CHECK(map.getNumOwnedElements("map1") == 3);
  CHECK(map.getNumOwnedElements("map2") == 5);

  // check num owned/used
  if (NumProc > 1) {
    if (MyPID == 0 || MyPID == NumProc - 1) {
      CHECK(map.getNumUsedElements("map1") == 4);
      CHECK(map.getNumUsedElements("map2") == 6);
    } else {
      CHECK(map.getNumUsedElements("map1") == 5);
      CHECK(map.getNumUsedElements("map2") == 7);
    }
  } else {
    CHECK(map.getNumUsedElements("map1") == 3);
    CHECK(map.getNumUsedElements("map2") == 5);
  }

  // check num dofs
  CHECK(map.getNumVectors("map1") == 2);
  CHECK(map.getNumVectors("map2") == 2);

  // check that the maps properly export
  Vector_type owned(map.getMap());
  Vector_type ghosted(map.getGhostedMap());
  Import_type importer(map.getMap(), map.getGhostedMap());

  for (int i = 0; i != map.getMap()->getLocalNumElements(); ++i) {
    owned.replaceLocalValue(i, map.getMap()->getGlobalElement(i));
  }

  ghosted.doImport(owned, importer, Tpetra::INSERT);
  ghosted.print(std::cout);
  {
    auto ghosted_view = ghosted.getData();
    for (int i = 0; i != map.getGhostedMap()->getLocalNumElements(); ++i) {
      CHECK(ghosted_view[i] == map.getGhostedMap()->getGlobalElement(i));
    }
  }

  // check the indices
  {
    auto inds_m1_d0 = map.viewIndices("map1", 0);
    CHECK(inds_m1_d0.size() == 3);
    CHECK(inds_m1_d0[0] == 0);
    CHECK(inds_m1_d0[1] == 2);
    CHECK(inds_m1_d0[2] == 4);

    auto inds_m1_d1 = map.viewIndices("map1", 1);
    CHECK(inds_m1_d1.size() == 3);
    CHECK(inds_m1_d1[0] == 1);
    CHECK(inds_m1_d1[1] == 3);
    CHECK(inds_m1_d1[2] == 5);

    auto inds_m2_d0 = map.viewIndices("map2", 0);
    CHECK(inds_m2_d0.size() == 5);
    CHECK(inds_m2_d0[0] == 6);
    CHECK(inds_m2_d0[1] == 8);
    CHECK(inds_m2_d0[2] == 10);
    CHECK(inds_m2_d0[3] == 12);
    CHECK(inds_m2_d0[4] == 14);

    auto inds_m2_d1 = map.viewIndices("map2", 1);
    CHECK(inds_m2_d1.size() == 5);
    CHECK(inds_m2_d1[0] == 7);
    CHECK(inds_m2_d1[1] == 9);
    CHECK(inds_m2_d1[2] == 11);
    CHECK(inds_m2_d1[3] == 13);
    CHECK(inds_m2_d1[4] == 15);
  }

  {
    auto inds_m1_d0 = map.viewGhostIndices("map1", 0);
    CHECK(inds_m1_d0[0] == 0);
    CHECK(inds_m1_d0[1] == 2);
    CHECK(inds_m1_d0[2] == 4);

    auto inds_m1_d1 = map.viewGhostIndices("map1", 1);
    CHECK(inds_m1_d1[0] == 1);
    CHECK(inds_m1_d1[1] == 3);
    CHECK(inds_m1_d1[2] == 5);

    auto inds_m2_d0 = map.viewGhostIndices("map2", 0);
    CHECK(inds_m2_d0[0] == 6);
    CHECK(inds_m2_d0[1] == 8);
    CHECK(inds_m2_d0[2] == 10);
    CHECK(inds_m2_d0[3] == 12);
    CHECK(inds_m2_d0[4] == 14);

    auto inds_m2_d1 = map.viewGhostIndices("map2", 1);
    CHECK(inds_m2_d1[0] == 7);
    CHECK(inds_m2_d1[1] == 9);
    CHECK(inds_m2_d1[2] == 11);
    CHECK(inds_m2_d1[3] == 13);
    CHECK(inds_m2_d1[4] == 15);

    if (NumProc > 1) {
      if (MyPID == 0 || MyPID == NumProc - 1) {
        CHECK(inds_m1_d0.size() == 4);
        CHECK(inds_m1_d1.size() == 4);
        CHECK(inds_m2_d0.size() == 6);
        CHECK(inds_m2_d1.size() == 6);

        CHECK(inds_m1_d0[3] == 16);
        CHECK(inds_m1_d1[3] == 17);
        CHECK(inds_m2_d0[5] == 18);
        CHECK(inds_m2_d1[5] == 19);
      } else {
        CHECK(inds_m1_d0.size() == 5);
        CHECK(inds_m1_d1.size() == 5);
        CHECK(inds_m2_d0.size() == 7);
        CHECK(inds_m2_d1.size() == 7);

        CHECK(inds_m1_d0[3] == 16);
        CHECK(inds_m1_d0[4] == 18);
        CHECK(inds_m1_d1[3] == 17);
        CHECK(inds_m1_d1[4] == 19);
        CHECK(inds_m2_d0[5] == 20);
        CHECK(inds_m2_d0[6] == 22);
        CHECK(inds_m2_d1[5] == 21);
        CHECK(inds_m2_d1[6] == 23);
      }
    } else {
      CHECK(inds_m1_d0.size() == 3);
      CHECK(inds_m1_d1.size() == 3);
      CHECK(inds_m2_d0.size() == 5);
      CHECK(inds_m2_d1.size() == 5);
    }
  }
}

TEST(SUPERMAP_MANUAL_TEST)
{
  SuperMap_Manual(false);
  SuperMap_Manual(true);
}

TEST(SUPERMAP_FROM_SINGLE_COMPOSITEVECTOR)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int MyPID = comm->getRank();
  if (MyPID == 0)
    std::cout
      << "Test: SuperMapLumped from 1 CompositeVector with multiple components and multiple dofs"
      << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cv;
  std::vector<std::string> names;
  names.push_back("cell");
  names.push_back("face");
  std::vector<int> dofs{ 2, 2 };
  std::vector<Entity_kind> locs;
  locs.push_back(CELL);
  locs.push_back(FACE);
  cv.SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a SuperMapLumped from this space
  Teuchos::RCP<SuperMap> map = createSuperMap(cv);

  int ncells_owned = mesh->getNumEntities(CELL, Parallel_kind::OWNED);
  int nfaces_owned = mesh->getNumEntities(FACE, Parallel_kind::OWNED);
  int ncells_used = mesh->getNumEntities(CELL, Parallel_kind::ALL);
  int nfaces_used = mesh->getNumEntities(FACE, Parallel_kind::ALL);

  // check basic sizes
  CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned, map->getMap()->getLocalNumElements());
  CHECK_EQUAL(2 * ncells_used + 2 * nfaces_used, map->getGhostedMap()->getLocalNumElements());

  // check CompMaps
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::CELL,false)->isSameAs(*map->getComponentMap(0, "cell")));
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::CELL,true)->isSameAs(*map->getComponentGhostedMap(0, "cell")));
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::FACE,false)->isSameAs(*map->getComponentMap(0, "face")));
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::FACE,true)->isSameAs(*map->getComponentGhostedMap(0, "face")));

  // check ordering is as expected
  const auto& inds_c0 = map->viewGhostIndices<MirrorHost>(0, "cell", 0);
  const auto& inds_c1 = map->viewGhostIndices<MirrorHost>(0, "cell", 1);
  const auto& inds_f0 = map->viewGhostIndices<MirrorHost>(0, "face", 0);
  const auto& inds_f1 = map->viewGhostIndices<MirrorHost>(0, "face", 1);

  const auto& inds_c0_owned = map->viewIndices<MirrorHost>(0, "cell", 0);
  const auto& inds_c1_owned = map->viewIndices<MirrorHost>(0, "cell", 1);
  const auto& inds_f0_owned = map->viewIndices<MirrorHost>(0, "face", 0);
  const auto& inds_f1_owned = map->viewIndices<MirrorHost>(0, "face", 1);

  // check owned list is a subset of the ghost list
  CHECK_OWNED_SUBSET_GHOST(inds_c0_owned, ncells_owned, inds_c0, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_c1_owned, ncells_owned, inds_c1, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_f0_owned, nfaces_owned, inds_f0, nfaces_used);
  CHECK_OWNED_SUBSET_GHOST(inds_f1_owned, nfaces_owned, inds_f1, nfaces_used);

  // check lists are non-overlapping and of sufficient size to cover all values
  std::vector<cVectorView_type_<MirrorHost,LO> const *> inds_owned = { &inds_c0_owned, &inds_c1_owned, &inds_f0_owned, &inds_f1_owned };
  std::vector<cVectorView_type_<MirrorHost,LO> const *> inds_ghost = { &inds_c0, &inds_c1, &inds_f0, &inds_f1 };
  CHECK_UNIQUE(inds_owned, 2 * ncells_owned + 2 * nfaces_owned);
  CHECK_UNIQUE(inds_ghost, 2 * ncells_used + 2 * nfaces_used);

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
    CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned + 2, inds_c0[ncells_owned + 1]);
    CHECK_EQUAL(2 * ncells_owned + 2 * nfaces_owned + 3, inds_c1[ncells_owned + 1]);

    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned, inds_f0[nfaces_owned]);
    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned + 1, inds_f1[nfaces_owned]);
    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned + 2, inds_f0[nfaces_owned + 1]);
    CHECK_EQUAL(2 * ncells_used + 2 * nfaces_owned + 3, inds_f1[nfaces_owned + 1]);
  }
}


TEST(SUPERMAP_FROM_SINGLE_COMPOSITEVECTOR_REPEATED_MAPS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int MyPID = comm->getRank();
  if (MyPID == 0)
    std::cout
      << "Test: SuperMapLumped from 1 CompositeVector with multiple components which include "
         "repeated maps, should still interleave as if this were one map and two dofs"
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
  Teuchos::RCP<SuperMap> map = createSuperMap(cv);

  int ncells_owned = mesh->getNumEntities(CELL, Parallel_kind::OWNED);
  int ncells_used = mesh->getNumEntities(CELL, Parallel_kind::ALL);

  // check basic sizes
  CHECK_EQUAL(2 * ncells_owned, map->getMap()->getLocalNumElements());
  CHECK_EQUAL(2 * ncells_used, map->getGhostedMap()->getLocalNumElements());

  // check CompMaps
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::CELL,false)->isSameAs(*map->getComponentMap(0, "cellA")));
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::CELL,true)->isSameAs(*map->getComponentGhostedMap(0, "cellA")));
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::CELL,false)->isSameAs(*map->getComponentMap(0, "cellB")));
  CHECK(mesh->getMap(AmanziMesh::Entity_kind::CELL,true)->isSameAs(*map->getComponentGhostedMap(0, "cellB")));

  // check ordering is as expected
  const auto& inds_c0 = map->viewGhostIndices(0, "cellA", 0);
  const auto& inds_c1 = map->viewGhostIndices(0, "cellB", 0);

  const auto& inds_c0_owned = map->viewIndices(0, "cellA", 0);
  const auto& inds_c1_owned = map->viewIndices(0, "cellB", 0);

  // check owned list is a subset of the ghost list
  CHECK_OWNED_SUBSET_GHOST(inds_c0_owned, ncells_owned, inds_c0, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_c1_owned, ncells_owned, inds_c1, ncells_used);

  // check lists are non-overlapping and of sufficient size to cover all values
  std::vector<cVectorView_type_<MirrorHost,LO> const *> inds_owned = { &inds_c0_owned, &inds_c1_owned };
  CHECK_UNIQUE(inds_owned, 2 * ncells_owned);
  std::vector<cVectorView_type_<MirrorHost,LO> const *> inds_ghost = { &inds_c0, &inds_c1 };
  CHECK_UNIQUE(inds_ghost, 2 * ncells_used);

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
  Teuchos::RCP<SuperMap> map2 = createSuperMap(cv2);

  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));
}


TEST(SUPERMAP_FROM_TWO_IDENTICAL_COMPOSITEVECTORS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int MyPID = comm->getRank();
  if (MyPID == 0)
    std::cout << "Test: SuperMapLumped from 2 CompositeVectors with same map, single dof is same "
                 "as 1 CompositeVector with two dofs"
              << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 1 });
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 1 });

  // create a SuperMapLumped from this space
  Teuchos::RCP<SuperMap> map = Teuchos::rcp(new SuperMap(comm, { cvA.CreateSpace().ptr(),cvB.CreateSpace().ptr() }));

  // now create another with one CV, 2 dofs
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 2 });
  Teuchos::RCP<SuperMap> map2 = createSuperMap(cv2);

  // same map!
  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));

  // same indices!
  CHECK(CHECK_EQUAL_VIEW(map->viewIndices<Amanzi::DefaultHost>(0, "cell", 0),
                         map2->viewIndices<Amanzi::DefaultHost>(0, "cell", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 0),
                         map2->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewIndices<Amanzi::DefaultHost>(1, "cell", 0),
                         map2->viewIndices<Amanzi::DefaultHost>(0, "cell", 1)));
  CHECK(CHECK_EQUAL_VIEW(map->viewGhostIndices<Amanzi::DefaultHost>(1, "cell", 0),
                         map2->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 1)));
}


TEST(SUPERMAP_FROM_CELL_PLUS_FACE_IS_CELLFACE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int MyPID = comm->getRank();
  if (MyPID == 0)
    std::cout << "Test: SuperMapLumped from 2 CompositeVectors with different maps, single dof is "
                 "same as 1 CompositeVector with two components"
              << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell" }, std::vector<Entity_kind>{ CELL }, std::vector<int>{ 1 });
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "face" }, std::vector<Entity_kind>{ FACE }, std::vector<int>{ 1 });
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents(std::vector<std::string>{ "cell", "face" },
                                                 std::vector<Entity_kind>{ CELL, FACE },
                                                 std::vector<int>{ 1, 1 });

  // create a SuperMapLumped from this space
  auto map = Teuchos::rcp(new SuperMap(comm, { cvA.CreateSpace().ptr(),cvB.CreateSpace().ptr() }));
  auto map2 = Teuchos::rcp(new SuperMap(comm, { cv2.CreateSpace().ptr() }));

  // same map!
  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));

  // same indices!
  CHECK(CHECK_EQUAL_VIEW(map->viewIndices<Amanzi::DefaultHost>(0, "cell", 0),
                         map2->viewIndices<Amanzi::DefaultHost>(0, "cell", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 0),
                         map2->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewIndices<Amanzi::DefaultHost>(1, "face", 0),
                         map2->viewIndices<Amanzi::DefaultHost>(0, "face", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewGhostIndices<Amanzi::DefaultHost>(1, "face", 0),
                         map2->viewGhostIndices<Amanzi::DefaultHost>(0, "face", 0)));
}


TEST(SUPERMAP_FROM_SAME_NAME_DIFFERENT_MAP)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int MyPID = comm->getRank();
  if (MyPID == 0)
    std::cout << "Test: SuperMapLumped from 2 CompositeVectors with different maps but the same "
                 "name, is same as 1 CompositeVector with two components"
              << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell" }, std::vector<Entity_kind>{ CELL }, std::vector<int>{ 1 });
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents(
    std::vector<std::string>{ "cell" }, std::vector<Entity_kind>{ FACE }, std::vector<int>{ 1 });
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents(std::vector<std::string>{ "cell", "face" },
                                                 std::vector<Entity_kind>{ CELL, FACE },
                                                 std::vector<int>{ 1, 1 });

  // create a SuperMapLumped from this space
  auto map = Teuchos::rcp(new SuperMap(comm, { cvA.CreateSpace().ptr(),cvB.CreateSpace().ptr() }));
  auto map2 = Teuchos::rcp(new SuperMap(comm, { cv2.CreateSpace().ptr() }));

  // same map!
  CHECK(map->getMap()->isSameAs(*map2->getMap()));
  CHECK(map->getGhostedMap()->isSameAs(*map2->getGhostedMap()));

  // same indices!
  CHECK(CHECK_EQUAL_VIEW(map->viewIndices<Amanzi::DefaultHost>(0, "cell", 0),
                         map2->viewIndices<Amanzi::DefaultHost>(0, "cell", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 0),
                         map2->viewGhostIndices<Amanzi::DefaultHost>(0, "cell", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewIndices<Amanzi::DefaultHost>(1, "cell", 0),
                         map2->viewIndices<Amanzi::DefaultHost>(0, "face", 0)));
  CHECK(CHECK_EQUAL_VIEW(map->viewGhostIndices<Amanzi::DefaultHost>(1, "cell", 0),
                         map2->viewGhostIndices<Amanzi::DefaultHost>(0, "face", 0)));
}

//
// TPETRA DOES NOT HAVE BLOCKMAPS!
//

// TEST(SUPERMAP_FROM_SAME_NAME_SAME_MAP_DIFFERENT_ELEMENTSIZE)
// {
//   //
//   // This test should be good for pressure in matrix + fracture?
//   //
//   using namespace Amanzi;
//   using namespace Amanzi::AmanziMesh;
//   using namespace Amanzi::AmanziGeometry;

//   auto comm = getDefaultComm();
//   int MyPID = comm->getRank();
//   if (MyPID == 0)
//     std::cout << "Test: SuperMapLumped from 2 CompositeVectors with same elements but different "
//                  "element sizes in the same compname"
//               << std::endl;

//   auto mesh = getMesh(comm);
//   int ncells = mesh->getNumEntities(CELL, Parallel_kind::OWNED);

//   // create a CVSpace
//   // -- we need to construct BlockMaps that include the variable block size
//   // -- first set up a vector of the block sizes
//   const auto& cell_map = mesh->getMap(AmanziMesh::Entity_kind::CELL,false);
//   IntVector_type block_size(cell_map);
//   block_size.putScalar(1);
//   block_size.replaceLocalValue(1, 2);

//   // -- create the owned point map and thereby the owned block map
//   auto total_block_size = 1 + cell_map->getLocalNumElements();
//   auto point_map = Teuchos::rcp(new Map_type(-1, total_block_size, 0, comm));
//   auto block_cell_map = Teuchos::rcp(new BlockMap_type(point_map, cell_map->getMyGlobalIndices(), block_size.getData()));

//   // -- now the ghosted block map -- first get block sizes on ghosted entities by importing
//   const auto& cell_mapg = mesh->getMap(AmanziMesh::Entity_kind::CELL,true);
//   IntVector_type block_sizeg(cell_mapg);
//   block_sizeg.doImport(block_size, mesh->getImporter(Entity_kind::CELL), Tpetra::INSERT);

//   // -- then get the first point of every block by importing
//   Vector_type_<Entity_GID> first_point_gid(cell_map);
//   Entity_ID point_id = 0;
//   {
//     auto block_size_d = block_size.getData();
//     for (int i=0; i!=cell_map.getLocalNumElements(); ++i) {
//       first_point_gid.replaceLocalValue(i, point_map->getGlobalElement(point_id));
//       point_id += block_size_d[0];
//     }
//   }
//   Vector_type_<Entity_GID> first_point_gid_ghosted(cell_mapg);
//   first_point_gid_ghosted.doImport(first_point_gid, importer, Tpetra::INSERT);

//   // -- now we can construct the ghosted block map
//   auto block_cell_map_g = Teuchos::rcp(new Epetra_BlockMap(
//             cell_mapg->getNumGlobalElements(),
//             cell_mapg->getMyGlobalIndices(),
//             first_point_gid_ghosted.getData(),
//             block_sizeg.getData(), 0, comm));

//   CompositeVectorSpace cvA;
//   cvA.SetMesh(mesh)->SetGhosted()->SetComponents({ "cell" }, { CELL }, { 1 });
//   CompositeVectorSpace cvB;
//   cvB.SetMesh(mesh)->SetGhosted()->SetComponents(
//     { "cell" }, { CELL }, { { "cell", block_cell_map } }, { { "cell", block_cell_map_g } }, { 1 });

//   // create a SuperMapLumped from this space
//   Teuchos::RCP<SuperMap> map = Teuchos::rcp(new SuperMap({ cvA, cvB }));

//   // not the same map!
//   CHECK_EQUAL(ncells, map->viewIndices(0, "cell", 0).size());
//   CHECK_EQUAL(block_cell_map->getLocalNumPoints(), map->viewIndices(1, "cell", 0).size());
//   CHECK_EQUAL(1 + ncells, block_cell_map->getLocalNumPoints());

//   {
//     const auto& inds1 = map->viewIndices(0, "cell", 0);
//     const auto& inds2 = map->viewIndices(1, "cell", 0);
//     CHECK_EQUAL(0, inds1[0]);
//     CHECK_EQUAL(1, inds1[1]);

//     CHECK_EQUAL(ncells, inds2[0]);
//     CHECK_EQUAL(ncells + 1, inds2[1]);
//   }
//   CHECK_UNIQUE(std::vector<const std::vector<int>*>{ &inds1, &inds2 }, 2 * ncells + 1);
// }


// TEST(SUPERMAP_FROM_TREEVECTOR)
// {
//   using namespace Amanzi;
//   using namespace Amanzi::AmanziMesh;
//   using namespace Amanzi::AmanziGeometry;

//   auto comm = getDefaultComm();
//   int MyPID = comm->getRank();

//   if (MyPID == 0) std::cout << "Test: FD like matrix, null off-proc assembly" << std::endl;

//   // read parameter list
//   std::string xmlFileName = "test/operator_convergence.xml";
//   Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
//   Teuchos::ParameterList plist = xmlreader.getParameters();

//   Amanzi::VerboseObject::global_hide_line_prefix = true;

//   // create a mesh
//   Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
//   Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

//   Preference pref;
//   pref.clear();
//   pref.push_back(Framework::MSTK);

//   MeshFactory meshfactory(comm, gm);
//   meshfactory.set_preference(pref);
//   Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
//   //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

//   // create a CVSpace
//   Teuchos::RCP<CompositeVectorSpace> cv = Teuchos::rcp(new CompositeVectorSpace());
//   std::vector<std::string> names;
//   names.push_back("cell");
//   names.push_back("face");
//   std::vector<int> dofs(2, 2);
//   std::vector<Entity_kind> locs;
//   locs.push_back(CELL);
//   locs.push_back(FACE);
//   cv->SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

//   // create a TV
//   TreeVectorSpace tv;
//   Teuchos::RCP<TreeVectorSpace> tv_ss1 = Teuchos::rcp(new TreeVectorSpace());
//   tv_ss1->SetData(cv);
//   Teuchos::RCP<TreeVectorSpace> tv_ss2 = Teuchos::rcp(new TreeVectorSpace());
//   tv_ss2->SetData(cv);
//   tv.PushBack(tv_ss1);
//   tv.PushBack(tv_ss2);

//   // create a SuperMapLumped from a singleton space
//   Teuchos::RCP<SuperMap> map_singleton = createSuperMap(*tv_ss1);
//   // create a SuperMapLumped from a tree space
//   Teuchos::RCP<SuperMap> map = createSuperMap(tv);
// }
