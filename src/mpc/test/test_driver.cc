#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "State_Old.hh"
#include "MPC.hh"

#include "MeshFactory.hh"
#include "Mesh.hh"



TEST(DRIVER) {


  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  std::string    xmlInFileName = "test/driver.xml";


  
  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;
  Teuchos::ParameterXMLFileReader xmlreader(xmlInFileName);
  driver_parameter_list = xmlreader.getParameters();
  


  // For now create one geometric model from all the regions in the spec

  Teuchos::ParameterList reg_params = driver_parameter_list.sublist("Regions");

  int spdim = 3;
  Amanzi::AmanziGeometry::GeometricModelPtr 
      geom_model_ptr( new Amanzi::AmanziGeometry::GeometricModel(spdim, reg_params, comm) );


  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = driver_parameter_list.sublist("Mesh");

  std::string mesh_class = mesh_parameter_list.get<string>("Mesh Class");

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  
  cout << "Using mesh framework " << mesh_class << endl;

  if (mesh_class == "MSTK")  
    {      
      Teuchos::ParameterList mstk_mesh_parameter_list = 
      	mesh_parameter_list.sublist("MSTK Mesh Parameters");
      
      string filename = mstk_mesh_parameter_list.get<string>("Exodus file name");

      Amanzi::AmanziMesh::FrameworkPreference pref;
      pref.clear();
      pref.push_back(Amanzi::AmanziMesh::MSTK);
      
      Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
      meshfactory.preference(pref);
 
      mesh = meshfactory(filename.c_str(), geom_model_ptr);      
    }
  else
    {
      throw std::exception();
    }


  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  // create the MPC
  Amanzi::MPC mpc(driver_parameter_list, mesh, comm, obs_data);
  
  mpc.cycle_driver();
  
  
  delete comm;
      
}


