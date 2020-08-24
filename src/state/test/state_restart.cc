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

  TEST(RESTART_DUMP_REQUIRED_INTERVAL) {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base","restartdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 3;
    csps[2] = 10;
    plist.set<Teuchos::Array<int> >("cycles start period stop", csps);

    auto comm = Amanzi::getDefaultComm();
    Amanzi::Checkpoint R(plist, comm);
    
    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    int cycles_[31] = {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                       0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) {
      cycles[ic] = R.DumpRequested(ic);
    }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }

  TEST(RESTART_DUMP_REQUIRED_INTERVAL_OPENENDED1) {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);    

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 3;
    csps[2] = -1;
    plist.set<Teuchos::Array<int> >("cycles start period stop", csps);    

    auto comm = Amanzi::getDefaultComm();
    Amanzi::Checkpoint R(plist, comm);

    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    int cycles_[31] = {1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
                       1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1};
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) {
      cycles[ic] = R.DumpRequested(ic);
    }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_INTERVAL_OPENENDED2) {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base","restartdump");
    plist.set<int>("file name digits", 5);        

    Teuchos::Array<int> csps(3);
    csps[0] = 5;
    csps[1] = 3;
    csps[2] = -1;
    plist.set<Teuchos::Array<int> >("cycles start period stop", csps);    

    auto comm = Amanzi::getDefaultComm();
    Amanzi::Checkpoint R(plist, comm);

    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    int cycles_[31] = {0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 
                       0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) {
      cycles[ic] = R.DumpRequested(ic);
    }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_CYCLES) {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);    

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 3;
    csps[2] = 10;
    plist.set<Teuchos::Array<int> >("cycles start period stop", csps);    
    
    Teuchos::Array<int> cyc(3);
    cyc[0] = 2;
    cyc[1] = 3;
    cyc[2] = 4;
    plist.set<Teuchos::Array<int> >("cycles",cyc);
    
    auto comm = Amanzi::getDefaultComm();
    Amanzi::Checkpoint R(plist, comm);

    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    int cycles_[31] = {1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
                       0, 0, 0, 0, 
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) {
      cycles[ic] = R.DumpRequested(ic);
    }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_TIMES) {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base","restartdump");
    plist.set<int>("file name digits", 5);    

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 3.0;
    tsps[2] = 10.0;
    plist.set<Teuchos::Array<double> >("times start period stop", tsps);    
    
    Teuchos::Array<double> tim(3);
    tim[0] = 2.0;
    tim[1] = 3.0;
    tim[2] = 4.0;
    plist.set<Teuchos::Array<double> >("times",tim);
    
    auto comm = Amanzi::getDefaultComm();
    Amanzi::Checkpoint R(plist, comm);

    // test the cycle stuff, the expected result is in cycles_ and 
    // we store the computed result in cycles
    int times_[31] = { 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
                        0, 0, 0, 0, 
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };
    int times[31];
    for (int ic = 0; ic<=30; ic+=1) {
      times[ic] = R.DumpRequested((double)ic);
    }
    CHECK_ARRAY_EQUAL(times_, times, 31);
  }


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
    
    S0->set_time(1.02);

    Teuchos::ParameterList checkpoint_list = plist.get<Teuchos::ParameterList>("checkpoint");
    auto R = Teuchos::rcp( new Amanzi::Checkpoint(checkpoint_list, comm));

    WriteCheckpoint(*R, *S0, 0.0);

    auto S1 = Teuchos::rcp(new Amanzi::State(state_list) );
    S1->RequireField("celldata")->SetMesh(Mesh)->SetGhosted(false)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);

    S1->Setup();
    S1->InitializeFields();    

    // fill with random data before reading checkpoint
    Epetra_MultiVector& s1p = *S1->GetFieldData("celldata", "state")->ViewComponent("cell", false); 
    s1p.Random();

    ReadCheckpoint(comm, *S1, "restartdump00000.h5");

    Epetra_MultiVector& s0p = *S0->GetFieldData("celldata", "state")->ViewComponent("cell", false);     

    // and compare with the original
    CHECK_EQUAL(S0->time(), S1->time());

    CHECK_EQUAL(s0p.MyLength(), s1p.MyLength()); 

    for (int i = 0; i < s1p.MyLength(); ++i) {
      CHECK_EQUAL(s1p[0][i], s0p[0][i]);
    }
  }
}

