/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "UnitTest++.h"

#include "Checkpoint.hh"
#include "MeshFactory.hh"
#include "State.hh"


SUITE(RESTART) {
  TEST(DUMP_DATA) {
    Teuchos::ParameterList plist;
    // plist.set<std::string>("File Name Base","restart_dump");
    // plist.set<int>("File Name Digits",4);
    // Teuchos::ParameterList& i1_ = plist.sublist("Cycle Data");

    // i1_.set<int>("Start",0);
    // i1_.set<int>("End",10);
    // i1_.set<int>("Interval",1);

    auto comm = Amanzi::getDefaultComm();

    std::string xmlFileName = "test/state_restart.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    plist = xmlreader.getParameters();

    Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

    Amanzi::AmanziMesh::Preference pref;
    pref.clear();
    pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
    pref.push_back(Amanzi::AmanziMesh::Framework::STK);
    pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);

    Amanzi::AmanziMesh::MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(pref);
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> Mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

    Teuchos::ParameterList state_list = plist.get<Teuchos::ParameterList>("state");
    // now populate the parameter list...
    auto S0 = Teuchos::rcp(new Amanzi::State(state_list) );

    S0->RegisterDomainMesh(Mesh);

    S0->RequireField("celldata")->SetMesh(Mesh)->SetGhosted(false)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);

    S0->Setup();
    S0->InitializeFields();
    S0->set_cycle(0);
    S0->set_time(1.02);

    Teuchos::ParameterList checkpoint_list = plist.get<Teuchos::ParameterList>("checkpoint");
    auto R = Teuchos::rcp( new Amanzi::Checkpoint(checkpoint_list, *S0));
    R->Write(*S0, 0.0);

    auto S1 = Teuchos::rcp(new Amanzi::State(state_list) );
    S1->RequireField("celldata")->SetMesh(Mesh)->SetGhosted(false)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);

    S1->RegisterDomainMesh(Mesh);
    S1->Setup();
    S1->InitializeFields();

    // fill with random data before reading checkpoint
    Epetra_MultiVector& s1p = *S1->GetFieldData("celldata", "state")->ViewComponent("cell", false);
    s1p.Random();

    ReadCheckpoint(*S1, "restartdump00000.h5");

    Epetra_MultiVector& s0p = *S0->GetFieldData("celldata", "state")->ViewComponent("cell", false);

    // and compare with the original
    CHECK_EQUAL(S0->time(), S1->time());

    CHECK_EQUAL(s0p.MyLength(), s1p.MyLength());

    for (int i = 0; i < s1p.MyLength(); ++i) {
      CHECK_EQUAL(s1p[0][i], s0p[0][i]);
    }
  }
}

