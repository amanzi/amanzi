/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

// TPLs
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CompositeVector.hh"
#include "MeshFactory.hh"

// Amanzi::State
#include "IO.hh"
#include "State.hh"
#include "Visualization.hh"

SUITE(VISUALIZATION)
{
  TEST(VIZ_DUMP_REQUIRED)
  {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "visdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 4;
    csps[2] = 10;
    plist.set<Teuchos::Array<int>>("cycles start period stop", csps);

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 4.0;
    tsps[2] = 10.0;
    plist.set<Teuchos::Array<double>>("times start period stop", tsps);

    Teuchos::Array<double> times(2);
    times[0] = 1.0;
    times[1] = 3.0;
    plist.set<Teuchos::Array<double>>("times", times);

    Epetra_MpiComm comm(MPI_COMM_WORLD);
    Amanzi::Visualization V(plist);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[31] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) { cycles[ic] = V.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);

    // test the time sps stuff
    CHECK_EQUAL(true, V.DumpRequested(0.0));
    CHECK_EQUAL(true, V.DumpRequested(1.0));
    CHECK_EQUAL(true, V.DumpRequested(3.0));
    CHECK_EQUAL(true, V.DumpRequested(4.0));
    CHECK_EQUAL(true, V.DumpRequested(8.0));

    CHECK_EQUAL(false, V.DumpRequested(0.5));
    CHECK_EQUAL(false, V.DumpRequested(1.1));
    CHECK_EQUAL(false, V.DumpRequested(3.2));
    CHECK_EQUAL(false, V.DumpRequested(3.99));
    CHECK_EQUAL(false, V.DumpRequested(10.0));
  }

  TEST(DUMP_MESH_AND_DATA)
  {
    // here we just check that the code does not crash when
    // the mesh and data files are written
    Teuchos::ParameterList plist;

    auto comm = Amanzi::getDefaultComm();

    std::string xmlFileName = "test/state_vis.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    plist = xmlreader.getParameters();

    Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
    auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

    Amanzi::AmanziMesh::Preference pref;
    pref.clear();
    pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);

    Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
    meshfactory.set_preference(pref);
    auto Mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

    auto state_list = plist.sublist("state");
    Amanzi::State S0(state_list);

    S0.RegisterMesh("domain", Mesh);
    S0.Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>("celldata",
                                                                      Amanzi::Tags::DEFAULT)
      .SetMesh(Mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S0.Setup();
    S0.InitializeFields();

    S0.set_time(1.02);

    Teuchos::ParameterList visualization_list = plist.get<Teuchos::ParameterList>("visualization");
    Amanzi::Visualization V(visualization_list);
    V.set_mesh(Mesh);
    V.CreateFiles();

    WriteVis(V, S0);
  }
}
