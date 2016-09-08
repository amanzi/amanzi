#include <iostream>
#include "stdlib.h"
#include "math.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "mpc_pks_registration.hh"
#include "pks_chemistry_registration.hh"
#include "pks_transport_registration.hh"
#include "State.hh"


void RunTestReactiveTransport(const std::string& xmlInFileName) {
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  // read the main parameter list
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, &comm));

  // create mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 100, 1, 1, gm);
  ASSERT(!mesh.is_null());

  // create dummy observation data object
  double avg1, avg2;
  Amanzi::ObservationData obs_data;    
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::rcp(new Teuchos::ParameterList(plist));
 
  {
    Amanzi::CycleDriver cycle_driver(glist, mesh, &comm, obs_data);
    try {
      auto S = cycle_driver.Go();
      S->GetFieldData("total_component_concentration")->MeanValue(&avg1);
    } catch (...) {
      CHECK(false);
    }
  }

  // restart simulation and compare results
  glist->sublist("cycle driver").sublist("restart").set<std::string>("file name", "chk_rt00005.h5");
  glist->sublist("state").sublist("initial conditions").remove("geochemical conditions", false);

  {
    Amanzi::CycleDriver cycle_driver(glist, mesh, &comm, obs_data);
    try {
      auto S = cycle_driver.Go();
      S->GetFieldData("total_component_concentration")->MeanValue(&avg2);
    } catch (...) {
      CHECK(false);
    }
  }

  CHECK_CLOSE(avg1, avg2, 1e-5 * avg1);
}


TEST(MPC_DRIVER_REACTIVE_TRANSPORT_NATIVE) {
  RunTestReactiveTransport("test/mpc_driver_reactive_transport.xml");
}

TEST(MPC_DRIVER_REACTIVE_TRANSPORT_ALQUIMIA) {
  RunTestReactiveTransport("test/mpc_driver_alquimia_transport.xml");
}
