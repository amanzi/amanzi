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
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "IO.hh"
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "State.hh"


TEST(MPC_DRIVER_FLOW)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();

  // read the main parameter list
  std::string xmlInFileName = "test/mpc_flow.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory.create(0.0, 0.0, 216.0, 120.0, 54, 30);
  AMANZI_ASSERT(!mesh.is_null());

  // create dummy observation data object
  double avg1, avg2;
  Amanzi::ObservationData obs_data;
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::rcp(new Teuchos::ParameterList(plist));

  Teuchos::ParameterList state_plist = glist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  {
    Amanzi::CycleDriver cycle_driver(glist, S, comm, obs_data);
    // try {
    cycle_driver.Go();
    S->Get<CompositeVector>("saturation_liquid").MeanValue(&avg1);
    // } catch (const std::exception& e) {
    //   std::cerr << e.what() << "\n\n";
    //   ;
    //   CHECK(false);
    // } catch (...) {
    //   CHECK(false);
    // }
  }
  WriteStateStatistics(*S);
  S = Teuchos::null;

  // restart simulation and compare results
  glist->sublist("cycle driver")
    .sublist("restart")
    .set<std::string>("file name", "chk_flow00010.h5");
  avg2 = 0.0;
  S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  {
    Amanzi::CycleDriver cycle_driver(glist, S, comm, obs_data);
    try {
      cycle_driver.Go();
      S->Get<CompositeVector>("saturation_liquid").MeanValue(&avg2);
    } catch (...) {
      CHECK(false);
    }
  }
  WriteStateStatistics(*S);

  std::cout << "avg = " << avg1 << "," << avg2 << std::endl;
  CHECK_CLOSE(avg1, avg2, 1e-6 * avg1);

  // checking that we created only two pks
  CHECK(PKFactory::num_pks == 2);
  std::cout << PKFactory::list_pks << std::endl;
}
