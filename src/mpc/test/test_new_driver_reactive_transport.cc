#include <iostream>
#include "stdlib.h"
#include "math.h"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "CycleDriver.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_reactivetransport_registration.hh"
#include "pks_transport_registration.hh"
#include "pks_chemistry_registration.hh"
#include "State.hh"


TEST(NEW_DRIVER_Reactive_Transport) {

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  std::string xmlInFileName = "test/test_new_driver_reactive_transport.xml";

  // read the main parameter list
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
  Teuchos::RCP<Mesh> mesh = meshfactory("test/rect2D_10x10_ss.exo", gm);
  ASSERT(!mesh.is_null());

  bool mpc_new = true;
  
  // if (input_parameter_list.isParameter("New multi-process coordinator")){
  //   mpc_new = input_parameter_list.get<bool>("New multi-process coordinator",false);
  //   //mpc_new = true;
  // }

  // create dummy observation data object
  Amanzi::ObservationData obs_data;    





   

  if (mpc_new){
    if (driver_parameter_list.isSublist("State")){
      // Create the state.    
      Teuchos::ParameterList state_plist = driver_parameter_list.sublist("State");
      Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
      S->RegisterMesh("domain",mesh);      

      // -------------- MULTI-PROCESS COORDINATOR------- --------------------
      Amanzi::CycleDriver cycle_driver(driver_parameter_list, S, comm, obs_data);
      //--------------- DO THE SIMULATION -----------------------------------
      cycle_driver.go();
      //-----------------------------------------------------
    }
  }

  
  delete comm;
}


