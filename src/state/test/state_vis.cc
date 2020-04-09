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

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "CompositeVector.hh"
#include "MeshFactory.hh"

// State
#include "State.hh"
#include "io/Visualization.hh"

SUITE(VISUALIZATION)
{
  TEST(DUMP_MESH_AND_DATA)
  {
    // here we just check that the code does not crash when
    // the mesh and data files are written
    Teuchos::ParameterList plist;
    auto comm = Amanzi::getDefaultComm();

    std::string xmlFileName = "test/state_vis.xml";
    Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
    plist = xmlreader.getParameters();

    Teuchos::ParameterList region_list =
      plist.get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
      new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

    Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
    auto mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

    auto state_list = plist.sublist("state");
    Amanzi::State S0(state_list);

    S0.RegisterMesh("domain", mesh);
    S0.Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>(
        "celldata")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S0.Setup();
    S0.InitializeFields();

    S0.set_time(1.02);

    Teuchos::ParameterList visualization_list =
      plist.get<Teuchos::ParameterList>("visualization");
    Amanzi::Visualization V(visualization_list, mesh);
    WriteVis(V, S0);
  }
}
