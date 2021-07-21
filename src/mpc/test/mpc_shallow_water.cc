#include <iostream>
#include "stdlib.h"
#include "math.h"


#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "CycleDriver.hh"
#include "eos_registration.hh"
#include "MeshFactory.hh"
#include "Mesh.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_shallow_water_registration.hh"
#include "State.hh"


TEST(MPC_DRIVER_SHALLOW_WATER) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  auto comm = Amanzi::getDefaultComm();
  
  // read the main parameter list
  std::string xmlInFileName = "test/mpc_shallow_water.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  // create mesh
  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory.create(0.0, 0.0, 10.0, 10.0, 20, 20);
  AMANZI_ASSERT(!mesh.is_null());

  // create dummy observation data object
  double vmin;
  Amanzi::ObservationData obs_data;    
  Teuchos::RCP<Teuchos::ParameterList> glist = Teuchos::rcp(new Teuchos::ParameterList(plist));

  Teuchos::ParameterList state_plist = glist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("surface", mesh);
  
  {
    Amanzi::CycleDriver cycle_driver(glist, S, comm, obs_data);
    try {
      cycle_driver.Go();
      S->GetFieldData("surface-ponded_depth")->MinValue(&vmin);
    } catch (...) {
      CHECK(false);
    }
  }
  S = Teuchos::null;

  CHECK(vmin > 0.0);  

  // checking that we created only one pk
  CHECK(PKFactory::num_pks == 1);
  std::cout << PKFactory::list_pks << std::endl;
}


