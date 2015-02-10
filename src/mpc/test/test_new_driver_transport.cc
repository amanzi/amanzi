#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "State.hh"
#include "CycleDriver.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_transport_registration.hh"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Domain.hh"
#include "GeometricModel.hh"


TEST(NEW_DRIVER_TRANSPORT) {

using namespace std;
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;



#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  std::string xmlInFileName = "test/test_new_driver_transport.xml";

  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  driver_parameter_list = xmlreader.getParameters();

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList reg_params = driver_parameter_list.sublist("Regions");

  int spdim = 2;
  Amanzi::AmanziGeometry::GeometricModelPtr 
      geom_model_ptr( new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, comm) );

  Amanzi::AmanziGeometry::Domain *simdomain_ptr = new Amanzi::AmanziGeometry::Domain(spdim);

  simdomain_ptr->Add_Geometric_Model(geom_model_ptr);

  // ---------------- MESH -----------------------------------------------
  int rank, ierr, aerr, size;

  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = driver_parameter_list.sublist("Mesh");


  Amanzi::VerboseObject *meshverbobj = 
    new Amanzi::VerboseObject("Mesh", driver_parameter_list);

  // Create a mesh factory for this geometric model
  Amanzi::AmanziMesh::MeshFactory factory(comm,meshverbobj) ;

  // get the Mesh sublist
  ierr = 0;
  Teuchos::ParameterList mesh_params = driver_parameter_list.sublist("Mesh");
  
  Teuchos::ParameterList unstr_mesh_params = mesh_params.sublist("Unstructured");

  // Decide on which mesh framework to use
  bool expert_params_specified = unstr_mesh_params.isSublist("Expert");

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
// Read or generate the mesh
  std::string file(""), format("");



  if (unstr_mesh_params.isSublist("Generate Mesh")) {  // If Read parameters are specified
    Teuchos::ParameterList gen_params = unstr_mesh_params.sublist("Generate Mesh");
    ierr = 0;
    
    try {
      // create the mesh by internal generation
      mesh = factory.create(gen_params, geom_model_ptr);

    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }
  
    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) {
      exit(-aerr);
    }

  }else if (unstr_mesh_params.isSublist("Read Mesh File")) {
    Teuchos::ParameterList read_params = unstr_mesh_params.sublist("Read Mesh File");
    
    if (read_params.isParameter("File")) {
      file = read_params.get<std::string>("File");
    } else {
      std::cerr << "Must specify File parameter for Read option under Mesh" << std::endl;
      throw std::exception();
    }

    if (read_params.isParameter("Format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("Format");

      if (format != "Exodus II" && format != "H5M") {	    
	std::cerr << "Can only read files in Exodus II or H5M format" << std::endl;
	throw std::exception();
      }
    } else {
      std::cerr << "Must specify Format parameter for Read option under Mesh" << std::endl;
      throw std::exception();
    }

    if (!file.empty()) {
      ierr = 0;
      try {
        // create the mesh from the file
	mesh = factory.create(file, geom_model_ptr);
	    
      } catch (const std::exception& e) {
	std::cerr << rank << ": error: " << e.what() << std::endl;
	ierr++;
      }
	  
      comm->SumAll(&ierr, &aerr, 1);
      if (aerr > 0) {
	exit(-aerr);
      }
    }

  } 
  else {  // If Generate parameters are specified
    std::cerr << rank << ": error: " << "Neither Read nor Generate options specified for mesh" << std::endl;
    throw std::exception();
  }

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


