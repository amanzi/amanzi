/*
  This is the operator component of the Amanzi code. 

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

#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"

#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"
#include "SuperMap.hh"

/* *****************************************************************
 * manually constructed test
 * **************************************************************** */
TEST(SUPERMAP_MANUAL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  //  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();
  int NumProc = comm->NumProc();

  if (MyPID == 0) std::cout << "Test: Manual test of SuperMap" << std::endl;

  // make a ghosted and local map 1
  std::vector<int> gids(3);
  for (int i=0; i!=3; ++i) {
    gids[i] = 3*MyPID + i;
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
  std::vector<Teuchos::RCP<const Epetra_Map> > maps; maps.push_back(owned_map1); maps.push_back(owned_map2);
  std::vector<Teuchos::RCP<const Epetra_Map> > gmaps; gmaps.push_back(ghosted_map1); gmaps.push_back(ghosted_map2);
  Operators::SuperMap map(comm, names, dofnums, maps, gmaps);

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
  
  
TEST(SUPERMAP_FROM_COMPOSITEVECTOR) {
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
  CompositeVectorSpace cv;
  std::vector<std::string> names; names.push_back("cell"); names.push_back("face");
  std::vector<int> dofs(2,2);
  std::vector<Entity_kind> locs; locs.push_back(CELL); locs.push_back(FACE);
  cv.SetMesh(mesh)->SetGhosted()->SetComponents(names, locs, dofs);

  // create a SuperMap from this space
  Teuchos::RCP<SuperMap> map = createSuperMap(cv);
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
  
  // create a SuperMap from a singleton space
  Teuchos::RCP<SuperMap> map_singleton = createSuperMap(*tv_ss1);
  // create a SuperMap from a tree space
  Teuchos::RCP<SuperMap> map = createSuperMap(tv);
}
