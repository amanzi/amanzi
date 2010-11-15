#include <iostream>
#include "stdlib.h"
#include "math.h"

#include "Mesh_maps_moab.hh"
#include "Mesh_maps_simple.hh"
#include "Mesh_maps_base.hh"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Version.hpp"


#include "State.hpp"
#include "MPC.hpp"


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  
  Teuchos::CommandLineProcessor CLP;
  
  CLP.setDocString("\nThe Amanzi driver reads an XML input file and\n"
		   "runs a reactive flow and transport simulation.\n");
  
  std::string xmlInFileName = "options.xml";
  CLP.setOption("xml_file", &xmlInFileName, "XML options file");
  
  CLP.throwExceptions(false);
  
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn = CLP.parse(argc, argv);
  
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return 1;
  }
  
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  
  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;
  Teuchos::updateParametersFromXmlFile(xmlInFileName,&driver_parameter_list);
  
  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = driver_parameter_list.sublist("Mesh");

  std::string mesh_class = mesh_parameter_list.get<string>("Mesh Class");

  Teuchos::RCP<Mesh_maps_base> mesh;
  
  if (mesh_class == "Simple") 
    {
      Teuchos::ParameterList simple_mesh_parameter_list = 
      	mesh_parameter_list.sublist("Simple Mesh Parameters");

      Teuchos::RCP<Mesh_maps_simple> MMS = 
      	Teuchos::rcp(new Mesh_maps_simple(simple_mesh_parameter_list, comm));
      
      mesh = MMS;
      
    } 
  else if (mesh_class == "MOAB")  
    {
      
      Teuchos::ParameterList moab_mesh_parameter_list = 
      	mesh_parameter_list.sublist("MOAB Mesh Parameters");
      
      string filename = moab_mesh_parameter_list.get<string>("Exodus file name");

      Teuchos::RCP<Mesh_maps_moab> MMM = 
      	Teuchos::rcp(new Mesh_maps_moab(filename.c_str(), MPI_COMM_WORLD));      
      
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


