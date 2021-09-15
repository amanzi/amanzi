/*
  Data Structures

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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

#include "Epetra_BlockMap.h"
#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#include "MeshFactory.hh"

#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"

#define SUPERMAP_TESTING 1
#include "SuperMapLumped.hh"

#include "SuperMap.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;


Teuchos::RCP<Mesh> getMesh(const Comm_ptr_type& comm) {
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

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  return meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
}



void CHECK_OWNED_SUBSET_GHOST(const std::vector<int>& owned, int n_owned,
        const std::vector<int>& ghost, int n_ghost) {
  CHECK_EQUAL(n_owned, owned.size());
  CHECK_EQUAL(n_ghost, ghost.size());
  for (int i=0; i!=n_owned; ++i) {
    CHECK_EQUAL(owned[i], ghost[i]);
  }
}

void CHECK_UNIQUE(std::vector<const std::vector<int>* > index_lists, int total) {
  std::set<int> all;
  int expected_count = 0;
  for (auto index_list : index_lists) {
    expected_count += index_list->size();
    all.insert(index_list->begin(), index_list->end());
  }
  CHECK_EQUAL(expected_count, all.size());
  CHECK_EQUAL(total, expected_count);
}


/* *****************************************************************
 * manually constructed test
 * **************************************************************** */
void SuperMap_Manual(bool continuous)
{
  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  int NumProc = comm->NumProc();

  if (MyPID == 0) std::cout << "Test: Manual test of SuperMapLumped" << std::endl;

  // make a ghosted and local map 1
  std::vector<int> gids(3);
  if (continuous) {
    for (int i=0; i!=3; ++i) {
      gids[i] = 3*MyPID + i;
    }
  } else {
    for (int i=0; i!=3; ++i) {
      gids[i] = 3*(NumProc - MyPID) + i;
    }
  }
  Teuchos::RCP<Epetra_Map> owned_map1 = Teuchos::rcp(new Epetra_Map(3*NumProc, 3, &gids[0], 0, *comm));

  if (MyPID > 0) gids.push_back(3*MyPID-1);
  if (MyPID < NumProc-1) gids.push_back(3*(MyPID+1));
  Teuchos::RCP<Epetra_Map> ghosted_map1 = Teuchos::rcp(new Epetra_Map(-1, gids.size(), &gids[0], 0, *comm));

  // make a ghosted and local map 2
  std::vector<int> gids2(5);
  for (int i=0; i!=5; ++i) {
    gids2[i] = 5*MyPID + i;
  }
  Teuchos::RCP<Epetra_Map> owned_map2 = Teuchos::rcp(new Epetra_Map(5*NumProc, 5, &gids2[0], 0, *comm));

  if (MyPID > 0) gids2.push_back(5*MyPID-1);
  if (MyPID < NumProc-1) gids2.push_back(5*(MyPID+1));
  Teuchos::RCP<Epetra_Map> ghosted_map2 = Teuchos::rcp(new Epetra_Map(-1, gids2.size(), &gids2[0], 0, *comm));

  // make the supermap
  std::vector<std::string> names; names.push_back("map1"); names.push_back("map2");
  std::vector<int> dofnums(2,2);

  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  maps.push_back(owned_map1);
  maps.push_back(owned_map2);
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > gmaps;
  gmaps.push_back(ghosted_map1);
  gmaps.push_back(ghosted_map2);

  Operators::SuperMapLumped map(comm, names, dofnums, maps, gmaps);

  // check the offsets
  CHECK(map.Offset("map1") == 0);
  CHECK(map.Offset("map2") == 2*3);

  // check the ghosted offsets
  CHECK(map.GhostedOffset("map1") == (2*3 + 2*5));
  if (NumProc > 1) {
    if (MyPID == 0 || MyPID == NumProc-1) {
      CHECK(map.GhostedOffset("map2") == (2*3 + 2*5 + 2));
    } else {
      CHECK(map.GhostedOffset("map2") == (2*3 + 2*5 + 4));
    }
  } else {
    CHECK(map.GhostedOffset("map2") == (2*3 + 2*5));
  }

  // check num owned
  CHECK(map.NumOwnedElements("map1") == 3);
  CHECK(map.NumOwnedElements("map2") == 5);

  // check num owned/used
  if (NumProc > 1) {
    if (MyPID == 0 || MyPID == NumProc-1) {
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
  CHECK(map.NumDofs("map1") == 2);
  CHECK(map.NumDofs("map2") == 2);

  // check that the maps properly export
  Epetra_Vector owned(*map.Map());
  Epetra_Vector ghosted(*map.GhostedMap());
  Epetra_Import importer(*map.GhostedMap(), *map.Map());
  
  for (int i=0; i!=map.Map()->NumMyElements(); ++i) {
    owned[i] = map.Map()->GID(i);
  }

  int ierr = ghosted.Import(owned, importer, Insert);
  ghosted.Print(std::cout);
  CHECK(!ierr);
  for (int i=0; i!=map.GhostedMap()->NumMyElements(); ++i) {
    CHECK(ghosted[i] == map.GhostedMap()->GID(i));
  }
  
  // check the indices
  {
    const std::vector<int>& inds_m1_d0 = map.Indices("map1", 0);
    CHECK(inds_m1_d0.size() == 3);
    CHECK(inds_m1_d0[0] == 0);
    CHECK(inds_m1_d0[1] == 2);
    CHECK(inds_m1_d0[2] == 4);
    
    const std::vector<int>& inds_m1_d1 = map.Indices("map1", 1);
    CHECK(inds_m1_d1.size() == 3);
    CHECK(inds_m1_d1[0] == 1);
    CHECK(inds_m1_d1[1] == 3);
    CHECK(inds_m1_d1[2] == 5);

    const std::vector<int>& inds_m2_d0 = map.Indices("map2", 0);
    CHECK(inds_m2_d0.size() == 5);
    CHECK(inds_m2_d0[0] == 6);
    CHECK(inds_m2_d0[1] == 8);
    CHECK(inds_m2_d0[2] == 10);
    CHECK(inds_m2_d0[3] == 12);
    CHECK(inds_m2_d0[4] == 14);

    const std::vector<int>& inds_m2_d1 = map.Indices("map2", 1);
    CHECK(inds_m2_d1.size() == 5);
    CHECK(inds_m2_d1[0] == 7);
    CHECK(inds_m2_d1[1] == 9);
    CHECK(inds_m2_d1[2] == 11);
    CHECK(inds_m2_d1[3] == 13);
    CHECK(inds_m2_d1[4] == 15);
  }

  {
    const std::vector<int>& inds_m1_d0 = map.GhostIndices("map1", 0);
    CHECK(inds_m1_d0[0] == 0);
    CHECK(inds_m1_d0[1] == 2);
    CHECK(inds_m1_d0[2] == 4);

    const std::vector<int>& inds_m1_d1 = map.GhostIndices("map1", 1);
    CHECK(inds_m1_d1[0] == 1);
    CHECK(inds_m1_d1[1] == 3);
    CHECK(inds_m1_d1[2] == 5);

    const std::vector<int>& inds_m2_d0 = map.GhostIndices("map2", 0);
    CHECK(inds_m2_d0[0] == 6);
    CHECK(inds_m2_d0[1] == 8);
    CHECK(inds_m2_d0[2] == 10);
    CHECK(inds_m2_d0[3] == 12);
    CHECK(inds_m2_d0[4] == 14);
    
    const std::vector<int>& inds_m2_d1 = map.GhostIndices("map2", 1);
    CHECK(inds_m2_d1[0] == 7);
    CHECK(inds_m2_d1[1] == 9);
    CHECK(inds_m2_d1[2] == 11);
    CHECK(inds_m2_d1[3] == 13);
    CHECK(inds_m2_d1[4] == 15);

    if (NumProc > 1) {
      if (MyPID == 0 || MyPID == NumProc-1) {
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
  
TEST(SUPERMAP_MANUAL_TEST) {
  SuperMap_Manual(false);
  SuperMap_Manual(true);
}
  
TEST(SUPERMAP_FROM_SINGLE_COMPOSITEVECTOR) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: SuperMapLumped from 1 CompositeVector with multiple components and multiple dofs" << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cv;
  std::vector<std::string> names; names.push_back("cell"); names.push_back("face");
  std::vector<int> dofs{2,2};
  std::vector<Entity_kind> locs; locs.push_back(CELL); locs.push_back(FACE);
  cv.SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = createSuperMap(cv);

  int ncells_owned = mesh->num_entities(CELL, Parallel_type::OWNED);
  int nfaces_owned = mesh->num_entities(FACE, Parallel_type::OWNED);
  int ncells_used = mesh->num_entities(CELL, Parallel_type::ALL);
  int nfaces_used = mesh->num_entities(FACE, Parallel_type::ALL);
  
  // check basic sizes
  CHECK_EQUAL(2*ncells_owned + 2*nfaces_owned, map->Map()->NumMyElements());
  CHECK_EQUAL(2*ncells_used + 2*nfaces_used, map->GhostedMap()->NumMyElements());

  // check CompMaps
  CHECK(mesh->cell_map(false).SameAs(*map->ComponentMap(0, "cell")));
  CHECK(mesh->cell_map(true).SameAs(*map->ComponentGhostedMap(0, "cell")));
  CHECK(mesh->face_map(false).SameAs(*map->ComponentMap(0, "face")));
  CHECK(mesh->face_map(true).SameAs(*map->ComponentGhostedMap(0, "face")));
  
  // check ordering is as expected
  const auto& inds_c0 = map->GhostIndices(0, "cell", 0);
  const auto& inds_c1 = map->GhostIndices(0, "cell", 1);
  const auto& inds_f0 = map->GhostIndices(0, "face", 0);
  const auto& inds_f1 = map->GhostIndices(0, "face", 1);

  const auto& inds_c0_owned = map->Indices(0, "cell", 0);
  const auto& inds_c1_owned = map->Indices(0, "cell", 1);
  const auto& inds_f0_owned = map->Indices(0, "face", 0);
  const auto& inds_f1_owned = map->Indices(0, "face", 1);

  // check owned list is a subset of the ghost list
  CHECK_OWNED_SUBSET_GHOST(inds_c0_owned, ncells_owned, inds_c0, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_c1_owned, ncells_owned, inds_c1, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_f0_owned, nfaces_owned, inds_f0, nfaces_used);
  CHECK_OWNED_SUBSET_GHOST(inds_f1_owned, nfaces_owned, inds_f1, nfaces_used);

  // check lists are non-overlapping and of sufficient size to cover all values
  CHECK_UNIQUE(std::vector<const std::vector<int>*>{ &inds_c0_owned, &inds_c1_owned, &inds_f0_owned, &inds_f1_owned }, 2*ncells_owned + 2*nfaces_owned);
  CHECK_UNIQUE(std::vector<const std::vector<int>*>{ &inds_c0, &inds_c1, &inds_f0, &inds_f1 }, 2*ncells_used + 2*nfaces_used);

  // check owned maps interleave
  CHECK_EQUAL(0, inds_c0[0]);
  CHECK_EQUAL(1, inds_c1[0]);
  CHECK_EQUAL(2, inds_c0[1]);
  CHECK_EQUAL(3, inds_c1[1]);

  CHECK_EQUAL(2*ncells_owned, inds_f0[0]);
  CHECK_EQUAL(2*ncells_owned+1, inds_f1[0]);
  CHECK_EQUAL(2*ncells_owned+2, inds_f0[1]);
  CHECK_EQUAL(2*ncells_owned+3, inds_f1[1]);


  // check ghosts pick up at the end of owned and interleave
  if (comm->NumProc() > 1) {
    CHECK_EQUAL(2*ncells_owned + 2*nfaces_owned, inds_c0[ncells_owned]);
    CHECK_EQUAL(2*ncells_owned + 2*nfaces_owned+1, inds_c1[ncells_owned]);
    CHECK_EQUAL(2*ncells_owned + 2*nfaces_owned+2, inds_c0[ncells_owned+1]);
    CHECK_EQUAL(2*ncells_owned + 2*nfaces_owned+3, inds_c1[ncells_owned+1]);

    CHECK_EQUAL(2*ncells_used + 2*nfaces_owned, inds_f0[nfaces_owned]);
    CHECK_EQUAL(2*ncells_used + 2*nfaces_owned+1, inds_f1[nfaces_owned]);
    CHECK_EQUAL(2*ncells_used + 2*nfaces_owned+2, inds_f0[nfaces_owned+1]);
    CHECK_EQUAL(2*ncells_used + 2*nfaces_owned+3, inds_f1[nfaces_owned+1]);
  }
}



TEST(SUPERMAP_FROM_SINGLE_COMPOSITEVECTOR_REPEATED_MAPS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: SuperMapLumped from 1 CompositeVector with multiple components which include repeated maps, should still interleave as if this were one map and two dofs" << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cv;
  std::vector<std::string> names; names.push_back("cellA"); names.push_back("cellB");
  std::vector<int> dofs{1,1};
  std::vector<Entity_kind> locs; locs.push_back(CELL); locs.push_back(CELL);
  cv.SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = createSuperMap(cv);

  int ncells_owned = mesh->num_entities(CELL, Parallel_type::OWNED);
  int ncells_used = mesh->num_entities(CELL, Parallel_type::ALL);
  
  // check basic sizes
  CHECK_EQUAL(2*ncells_owned, map->Map()->NumMyElements());
  CHECK_EQUAL(2*ncells_used, map->GhostedMap()->NumMyElements());

  // check CompMaps
  CHECK(mesh->cell_map(false).SameAs(*map->ComponentMap(0, "cellA")));
  CHECK(mesh->cell_map(true).SameAs(*map->ComponentGhostedMap(0, "cellA")));
  CHECK(mesh->cell_map(false).SameAs(*map->ComponentMap(0, "cellB")));
  CHECK(mesh->cell_map(true).SameAs(*map->ComponentGhostedMap(0, "cellB")));
  
  // check ordering is as expected
  const auto& inds_c0 = map->GhostIndices(0, "cellA", 0);
  const auto& inds_c1 = map->GhostIndices(0, "cellB", 0);

  const auto& inds_c0_owned = map->Indices(0, "cellA", 0);
  const auto& inds_c1_owned = map->Indices(0, "cellB", 0);

  // check owned list is a subset of the ghost list
  CHECK_OWNED_SUBSET_GHOST(inds_c0_owned, ncells_owned, inds_c0, ncells_used);
  CHECK_OWNED_SUBSET_GHOST(inds_c1_owned, ncells_owned, inds_c1, ncells_used);

  // check lists are non-overlapping and of sufficient size to cover all values
  CHECK_UNIQUE(std::vector<const std::vector<int>*>{ &inds_c0_owned, &inds_c1_owned }, 2*ncells_owned);
  CHECK_UNIQUE(std::vector<const std::vector<int>*>{ &inds_c0, &inds_c1 }, 2*ncells_used);

  // check owned maps interleave
  CHECK_EQUAL(0, inds_c0[0]);
  CHECK_EQUAL(1, inds_c1[0]);
  CHECK_EQUAL(2, inds_c0[1]);
  CHECK_EQUAL(3, inds_c1[1]);

  // check ghosts pick up at the end of owned and interleave
  if (comm->NumProc() > 1) {
    CHECK_EQUAL(2*ncells_owned, inds_c0[ncells_owned]);
    CHECK_EQUAL(2*ncells_owned+1, inds_c1[ncells_owned]);
    CHECK_EQUAL(2*ncells_owned+2, inds_c0[ncells_owned+1]);
    CHECK_EQUAL(2*ncells_owned+3, inds_c1[ncells_owned+1]);
  }

  // this map should be IDENTICAL to that of single component, 2 dofs
  CompositeVectorSpace cv2;
  std::vector<std::string> names2{"cell"};
  std::vector<int> dofs2{2};
  std::vector<Entity_kind> locs2{CELL};
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents(names2, locs2, dofs2);

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map2 = createSuperMap(cv2);

  CHECK(map->Map()->SameAs(*map2->Map()));
  CHECK(map->GhostedMap()->SameAs(*map2->GhostedMap()));
}


TEST(SUPERMAP_FROM_TWO_IDENTICAL_COMPOSITEVECTORS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: SuperMapLumped from 2 CompositeVectors with same map, single dof is same as 1 CompositeVector with two dofs" << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL}, {1});
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL}, {1});

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(new SuperMap({cvA, cvB}));

  // now create another with one CV, 2 dofs
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL}, {2});
  Teuchos::RCP<Operators::SuperMap> map2 = createSuperMap(cv2);

  // same map!
  CHECK(map->Map()->SameAs(*map2->Map()));
  CHECK(map->GhostedMap()->SameAs(*map2->GhostedMap()));

  // same indices!
  CHECK(map->Indices(0, "cell", 0) == map2->Indices(0, "cell", 0));
  CHECK(map->GhostIndices(0, "cell", 0) == map2->GhostIndices(0, "cell", 0));
  CHECK(map->Indices(1, "cell", 0) == map2->Indices(0, "cell", 1));
  CHECK(map->GhostIndices(1, "cell", 0) == map2->GhostIndices(0, "cell", 1));
}  


TEST(SUPERMAP_FROM_CELL_PLUS_FACE_IS_CELLFACE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: SuperMapLumped from 2 CompositeVectors with different maps, single dof is same as 1 CompositeVector with two components" << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents( std::vector<std::string>{"cell"}, std::vector<Entity_kind>{CELL}, std::vector<int>{1} );
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents( std::vector<std::string>{"face"}, std::vector<Entity_kind>{FACE}, std::vector<int>{1} );
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents( std::vector<std::string>{"cell","face"}, std::vector<Entity_kind>{CELL,FACE}, std::vector<int>{1,1} );

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(new SuperMap({cvA, cvB}));
  Teuchos::RCP<Operators::SuperMap> map2 = Teuchos::rcp(new SuperMap({cv2}));

  // same map!
  CHECK(map->Map()->SameAs(*map2->Map()));
  CHECK(map->GhostedMap()->SameAs(*map2->GhostedMap()));

  // same indices!
  CHECK(map->Indices(0, "cell", 0) == map2->Indices(0, "cell", 0));
  CHECK(map->GhostIndices(0, "cell", 0) == map2->GhostIndices(0, "cell", 0));
  CHECK(map->Indices(1, "face", 0) == map2->Indices(0, "face", 0));
  CHECK(map->GhostIndices(1, "face", 0) == map2->GhostIndices(0, "face", 0));
}  


TEST(SUPERMAP_FROM_SAME_NAME_DIFFERENT_MAP) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: SuperMapLumped from 2 CompositeVectors with different maps but the same name, is same as 1 CompositeVector with two components" << std::endl;

  auto mesh = getMesh(comm);

  // create a CVSpace
  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents( std::vector<std::string>{"cell"}, std::vector<Entity_kind>{CELL}, std::vector<int>{1} );
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents( std::vector<std::string>{"cell"}, std::vector<Entity_kind>{FACE}, std::vector<int>{1} );
  CompositeVectorSpace cv2;
  cv2.SetMesh(mesh)->SetGhosted()->SetComponents( std::vector<std::string>{"cell","face"}, std::vector<Entity_kind>{CELL,FACE}, std::vector<int>{1,1} );

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(new SuperMap({cvA, cvB}));
  Teuchos::RCP<Operators::SuperMap> map2 = Teuchos::rcp(new SuperMap({cv2}));

  // same map!
  CHECK(map->Map()->SameAs(*map2->Map()));
  CHECK(map->GhostedMap()->SameAs(*map2->GhostedMap()));

  // same indices!
  CHECK(map->Indices(0, "cell", 0) == map2->Indices(0, "cell", 0));
  CHECK(map->GhostIndices(0, "cell", 0) == map2->GhostIndices(0, "cell", 0));
  CHECK(map->Indices(1, "cell", 0) == map2->Indices(0, "face", 0));
  CHECK(map->GhostIndices(1, "cell", 0) == map2->GhostIndices(0, "face", 0));
}  


TEST(SUPERMAP_FROM_SAME_NAME_SAME_MAP_DIFFERENT_ELEMENTSIZE) {
  //
  // This test should be good for pressure in matrix + fracture?
  //
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: SuperMapLumped from 2 CompositeVectors with same elements but different element sizes in the same compname" << std::endl;

  auto mesh = getMesh(comm);
  int ncells = mesh->num_entities(CELL, Parallel_type::OWNED);

  // create a CVSpace
  const auto& cell_map = mesh->cell_map(false);
  const int * gids = nullptr;
  const long long * llgids = nullptr;
  cell_map.MyGlobalElements(gids, llgids);
  Epetra_IntVector element_size(cell_map);
  element_size.PutValue(1);
  element_size[1] = 2;

  Teuchos::RCP<const Epetra_BlockMap> block_cell_map = Teuchos::rcp(new Epetra_BlockMap(cell_map.NumGlobalElements(), cell_map.NumMyElements(), gids, &element_size[0], 0, *comm));

  const auto& cell_mapg = mesh->cell_map(true);
  Epetra_Import importer(cell_mapg, cell_map);
  Epetra_IntVector element_sizeg(cell_mapg);
  element_sizeg.Import(element_size, importer, Insert);
  
  cell_mapg.MyGlobalElements(gids, llgids);
  Teuchos::RCP<const Epetra_BlockMap> block_cell_map_g = Teuchos::rcp(new Epetra_BlockMap(cell_mapg.NumGlobalElements(), cell_mapg.NumMyElements(), gids, &element_sizeg[0], 0, *comm));
  

  CompositeVectorSpace cvA;
  cvA.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL}, {1});
  CompositeVectorSpace cvB;
  cvB.SetMesh(mesh)->SetGhosted()->SetComponents({"cell"}, {CELL}, {{"cell",block_cell_map}}, {{"cell", block_cell_map_g}}, {1});

  // create a SuperMapLumped from this space
  Teuchos::RCP<Operators::SuperMap> map = Teuchos::rcp(new SuperMap({cvA, cvB}));

  // not the same map!
  CHECK(map->Indices(0, "cell", 0).size() == ncells);
  CHECK(map->Indices(1, "cell", 0).size() == block_cell_map->NumMyPoints());
  CHECK_EQUAL(block_cell_map->NumMyPoints(), 1 + ncells);

  const auto& inds1 = map->Indices(0, "cell", 0);
  const auto& inds2 = map->Indices(1, "cell", 0);
  CHECK_EQUAL(0, inds1[0]);
  CHECK_EQUAL(1, inds1[1]);

  CHECK_EQUAL(ncells, inds2[0]);
  CHECK_EQUAL(ncells+1, inds2[1]);

  CHECK_UNIQUE(std::vector<const std::vector<int>*>{ &inds1, &inds2 }, 2*ncells + 1);
  
}  


TEST(SUPERMAP_FROM_TREEVECTOR) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: FD like matrix, null off-proc assembly" << std::endl;

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

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // create a CVSpace
  Teuchos::RCP<CompositeVectorSpace> cv = Teuchos::rcp(new CompositeVectorSpace());
  std::vector<std::string> names; names.push_back("cell"); names.push_back("face");
  std::vector<int> dofs(2,2);
  std::vector<Entity_kind> locs; locs.push_back(CELL); locs.push_back(FACE);
  cv->SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a TV
  TreeVectorSpace tv;
  Teuchos::RCP<TreeVectorSpace> tv_ss1 = Teuchos::rcp(new TreeVectorSpace());
  tv_ss1->SetData(cv);
  Teuchos::RCP<TreeVectorSpace> tv_ss2 = Teuchos::rcp(new TreeVectorSpace());
  tv_ss2->SetData(cv);
  tv.PushBack(tv_ss1);
  tv.PushBack(tv_ss2);
  
  // create a SuperMapLumped from a singleton space
  Teuchos::RCP<Operators::SuperMap> map_singleton = createSuperMap(*tv_ss1);
  // create a SuperMapLumped from a tree space
  Teuchos::RCP<Operators::SuperMap> map = createSuperMap(tv);
}
