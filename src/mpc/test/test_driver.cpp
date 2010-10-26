#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "Mesh_maps_simple.hh"
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

  Teuchos::RCP<Mesh_maps_simple> MMS;
  if (mesh_class == "Simple") 
    {
      Teuchos::ParameterList simple_mesh_parameter_list = 
	mesh_parameter_list.sublist("Simple Mesh Parameters");

      MMS = Teuchos::rcp(new Mesh_maps_simple(simple_mesh_parameter_list, comm));
      
    } 
  else 
    {
      // here's where we'll put the stuff for reading
      // parameter lists for STK and MOAB
      
      // for now just throw and exception
      throw std::exception();
    }


  // create the MPC
  MPC mpc(driver_parameter_list, MMS);
  
  mpc.write_mesh();

  mpc.cycle_driver();
  
  

      
}


