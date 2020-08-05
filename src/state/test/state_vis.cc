/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

// TPLs
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"

// State
#include "State.hh"
#include "Visualization.hh"


SUITE(VISUALIZATION) {
  TEST(DUMP_MESH_AND_DATA) {
    auto comm = Amanzi::getDefaultComm();

    std::string xmlFileName = "test/state_vis.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();

    Teuchos::ParameterList region_list = plist.sublist("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

    Amanzi::AmanziMesh::Preference pref;
    pref.clear();
    pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);   
    
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(pref);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> Mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

    Teuchos::ParameterList state_list = plist.sublist("state");
    Teuchos::RCP<Amanzi::State> S0 = Teuchos::rcp(new Amanzi::State(state_list));

    S0->RegisterMesh("domain",Mesh);
    S0->RequireField("celldata")->SetMesh(Mesh)->SetGhosted(false)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S0->Setup();
    S0->InitializeFields();
    
    S0->set_time(1.02);

    Teuchos::ParameterList visualization_list = plist.sublist("visualization");
    auto V = Teuchos::rcp(new Amanzi::Visualization(visualization_list));
    V->set_mesh(Mesh);
    V->CreateFiles();
    WriteVis(*V, *S0);
  }
}

