/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Markus Berndt
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

//#include "Epetra_MpiComm.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "CompositeVector.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "io/Checkpoint.hh"

SUITE(RESTART)
{
  TEST(DUMP_DATA)
  {
    std::string xmlFileName = "test/state_restart.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    Teuchos::ParameterList plist = xmlreader.getParameters();

    auto comm = Amanzi::getDefaultComm();

    // set up mesh
    Teuchos::ParameterList region_list = plist.sublist("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
      new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

    Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
    auto mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

    // set up state
    Teuchos::ParameterList state_list = plist.sublist("state");
    Amanzi::State S0(state_list);

    S0.RegisterDomainMesh(mesh);

    S0.Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>(
        "celldata", "", "state_restart")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", Amanzi::AmanziMesh::CELL, 3);

    S0.Setup();
    // randomize
    S0.GetW<Amanzi::CompositeVector>("celldata", "", "state_restart").random();
    S0.set_time(1.02);
    S0.set_cycle(5);

    // write checkpoint
    Teuchos::ParameterList checkpoint_list = plist.sublist("checkpoint");
    Amanzi::Checkpoint ckp(checkpoint_list, comm);
    WriteCheckpoint(ckp, S0);

    // read checkpoint
    Amanzi::State S1(state_list);
    S1.Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>(
        "celldata", "", "state_restart")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", Amanzi::AmanziMesh::CELL, 3);
    S1.Setup();

    S1.GetW<Amanzi::CompositeVector>("celldata", "state_restart").putScalar(0.);

    // -- read
    ReadCheckpoint(comm, S1, "restartdump00005.h5");

    // Compare
    auto s0p = S0.Get<Amanzi::CompositeVector>("celldata")
                 .ViewComponent<AmanziDefaultHost>("cell", false);
    auto s1p = S1.Get<Amanzi::CompositeVector>("celldata")
                 .ViewComponent<AmanziDefaultHost>("cell", false);

    // and compare with the original
    CHECK_EQUAL(S0.time(), S1.time());
    CHECK_EQUAL(s0p.extent(0), s1p.extent(0));
    CHECK_EQUAL(s0p.extent(1), s1p.extent(1));

    for (int i = 0; i < s0p.extent(0); ++i) {
      for (int j = 0; j < s0p.extent(1); ++j) {
        CHECK_EQUAL(s0p(i, j), s1p(i, j));
      }
    }
  }
}
