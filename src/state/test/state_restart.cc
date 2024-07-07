/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  State

*/

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Checkpoint.hh"
#include "CompositeVector.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "Tag.hh"

using namespace Amanzi;

SUITE(RESTART)
{
  TEST(RESTART_DUMP_REQUIRED_INTERVAL)
  {
    Teuchos::ParameterList plist;
    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 3;
    csps[2] = 10;
    plist.set<Teuchos::Array<int>>("cycles start period stop", csps);

    State S;
    Checkpoint R(plist, S);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[31] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) { cycles[ic] = R.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_INTERVAL_OPENENDED1)
  {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 3;
    csps[2] = -1;
    plist.set<Teuchos::Array<int>>("cycles start period stop", csps);

    State S;
    Checkpoint R(plist, S);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[31] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
                        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1 };
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) { cycles[ic] = R.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_INTERVAL_OPENENDED2)
  {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<int> csps(3);
    csps[0] = 5;
    csps[1] = 3;
    csps[2] = -1;
    plist.set<Teuchos::Array<int>>("cycles start period stop", csps);

    State S;
    Checkpoint R(plist, S);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[31] = { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 };
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) { cycles[ic] = R.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_CYCLES)
  {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 3;
    csps[2] = 10;
    plist.set<Teuchos::Array<int>>("cycles start period stop", csps);

    Teuchos::Array<int> cyc(3);
    cyc[0] = 2;
    cyc[1] = 3;
    cyc[2] = 4;
    plist.set<Teuchos::Array<int>>("cycles", cyc);

    State S;
    Checkpoint R(plist, S);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[31] = { 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) { cycles[ic] = R.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(RESTART_DUMP_REQUIRED_TIMES)
  {
    Teuchos::ParameterList plist;

    plist.set<std::string>("file name base", "restartdump");
    plist.set<int>("file name digits", 5);

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 3.0;
    tsps[2] = 10.0;
    plist.set<Teuchos::Array<double>>("times start period stop", tsps);

    Teuchos::Array<double> tim(3);
    tim[0] = 2.0;
    tim[1] = 3.0;
    tim[2] = 4.0;
    plist.set<Teuchos::Array<double>>("times", tim);

    State S;
    Checkpoint R(plist, S);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int times_[31] = { 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int times[31];
    for (int ic = 0; ic <= 30; ic += 1) { times[ic] = R.DumpRequested((double)ic); }
    CHECK_ARRAY_EQUAL(times_, times, 31);
  }


  TEST(DUMP_DATA)
  {
    auto comm = getDefaultComm();

    std::string xmlFileName = "test/state_restart.xml";
    auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

    Teuchos::ParameterList& region_list = plist->sublist("regions");
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

    AmanziMesh::Preference pref;
    pref.clear();
    pref.push_back(AmanziMesh::Framework::MSTK);
    pref.push_back(AmanziMesh::Framework::SIMPLE);

    AmanziMesh::MeshFactory meshfactory(comm, gm);
    meshfactory.set_preference(pref);
    auto Mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

    auto state_list = Teuchos::sublist(plist, "state");
    State S0(state_list);

    S0.RegisterDomainMesh(Mesh);

    S0.Require<CompositeVector, CompositeVectorSpace>("celldata", Tags::DEFAULT, "state_restart")
      .SetMesh(Mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    S0.Setup();
    S0.set_time(1.02);

    Teuchos::ParameterList checkpoint_list = plist->get<Teuchos::ParameterList>("checkpoint");
    Checkpoint R(checkpoint_list, S0);

    R.write(S0);

    State S1(state_list);
    S1.RegisterDomainMesh(Mesh);

    S1.Require<CompositeVector, CompositeVectorSpace>("celldata", Tags::DEFAULT, "state_restart")
      .SetMesh(Mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    S1.Setup();

    // fill with random data before reading checkpoint
    auto& s1p = *S1.GetW<CompositeVector>("celldata", "state_restart").getComponent("cell", false);
    s1p.randomize();

    Checkpoint chkp("restartdump00000.h5", comm);
    chkp.read(S1);

    auto& s0p = *S0.GetW<CompositeVector>("celldata", "state_restart").getComponent("cell", false);

    // and compare with the original
    CHECK_EQUAL(S0.get_time(), S1.get_time());
    CHECK_EQUAL(s0p.getLocalLength(), s1p.getLocalLength());

    s0p.update(-1, s1p, 1.);
    CHECK_CLOSE(0, s0p.getVector(0)->normInf(), 1.e-10);
  }
}
