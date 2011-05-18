#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"

#include "Mesh_MOAB.hh"
#include "Mesh_simple.hh"
#include "Mesh.hh"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "State.hpp"
#include "MPC.hpp"


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
  Teuchos::updateParametersFromXmlFile(xmlInFileName,&driver_parameter_list);
  
  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = driver_parameter_list.sublist("Mesh");

  std::string mesh_class = mesh_parameter_list.get<string>("Mesh Class");

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  
  cout << mesh_class << endl;

  if (mesh_class == "Simple") 
    {
      Teuchos::ParameterList simple_mesh_parameter_list = 
      	mesh_parameter_list.sublist("Simple Mesh Parameters");

      Teuchos::RCP<Amanzi::AmanziMesh::Mesh_simple> MMS = 
      	Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_simple(simple_mesh_parameter_list, comm));
      
      mesh = MMS;
      
    } 
  else if (mesh_class == "MOAB")  
    {
      
      Teuchos::ParameterList moab_mesh_parameter_list = 
      	mesh_parameter_list.sublist("MOAB Mesh Parameters");
      
      string filename = moab_mesh_parameter_list.get<string>("Exodus file name");

      Teuchos::RCP<Amanzi::AmanziMesh::Mesh_MOAB> MMM = 
      	Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MOAB(filename.c_str(), MPI_COMM_WORLD));      
      
      mesh = MMM;

    }
  else
    {
      throw std::exception();
    }


  // create the MPC
  MPC mpc(driver_parameter_list, mesh);
  
  mpc.cycle_driver();
  
  
  delete comm;
      
}


