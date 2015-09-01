#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "CycleDriver.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_registration.hh"
#include "pks_transport_registration.hh"
#include "pks_chemistry_registration.hh"
#include "State.hh"


TEST(MPC_DRIVER_FLOW_REACTIVE_TRANSPORT) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
  
  // read the main parameter list
  std::string xmlInFileName = "test/mpc_driver_frt.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);
   
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
  Amanzi::ObservationData obs_data;    
  Teuchos::RCP<Teuchos::ParameterList> glist_rcp = Teuchos::rcp(new Teuchos::ParameterList(plist));

  Amanzi::CycleDriver cycle_driver(glist_rcp, mesh, &comm, obs_data);
  cycle_driver.Go();

  // if (plist.isSublist("State")) {
  //   Teuchos::ParameterList state_plist = plist.sublist("State");
  //   Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  //   S->RegisterMesh("domain", mesh);      

  //   Amanzi::CycleDriver cycle_driver(glist_rcp, S, &comm, obs_data);
  //   cycle_driver.Go();
  // }
}


