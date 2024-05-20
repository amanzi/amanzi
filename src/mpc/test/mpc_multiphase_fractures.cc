/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "stdlib.h"
#include "math.h"

#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "IO.hh"
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "evaluators_multiphase_reg.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "models_multiphase_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_multiphase_reg.hh"
#include "State.hh"


TEST(MPC_DRIVER_MULTIPHASE_FRACTURES)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  // read the main parameter list
  std::string xmlInFileName = "test/mpc_multiphase_fractures.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");

  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh3D = factory.create(0.0, 0.0, 0.0, 200.0, 12.0, 60.0, 100, 2, 20);

  std::vector<std::string> names;
  names.push_back("fractures");
  auto mesh = factory.create(mesh3D, names, AmanziMesh::FACE);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  // work-around
  Key key("mass_density_gas");
  S->Require<CompositeVector, CompositeVectorSpace>(key, Tags::DEFAULT, key)
    .SetMesh(mesh)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  {
    Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
    try {
      cycle_driver.Go();
    } catch (const std::exception& e) {
      std::cerr << e.what() << "\n\n";
      CHECK(false);
    } catch (...) {
      CHECK(false);
    }
  }
  WriteStateStatistics(*S);
}
