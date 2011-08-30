#include <iostream>
#include <cstdlib>
#include <cmath>

#if HAVE_MOAB_MESH
#include "Mesh_MOAB.hh"
#endif

#if HAVE_STK_MESH
#include "MeshFactory.hh"
#include "Mesh_STK.hh"
#endif

#include "Mesh_simple.hh"
#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"
#include "Mesh.hh"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_DefaultMpiComm.hpp"


#include "State.hpp"
#include "MPC.hpp"

#include "errors.hh"
#include "exceptions.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);

  // make sure only PE0 can write to std::cout
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank!=0) {
    cout.rdbuf(0);
  } 

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
  
  // print parameter list
  std::cout << "======================> dumping parameter list <======================" << std::endl;
  Teuchos::writeParameterListToXmlOStream(driver_parameter_list, std::cout);
  std::cout << "======================> done dumping parameter list. <================"<<std::endl;

  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = driver_parameter_list.sublist("Mesh");

  std::string mesh_class = mesh_parameter_list.get<string>("Mesh Class");

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  
  if (mesh_class == "Simple") 
    {
      Teuchos::ParameterList simple_mesh_parameter_list = 
      	mesh_parameter_list.sublist("Simple Mesh Parameters");

      mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_simple(simple_mesh_parameter_list, comm));
    } 
#ifdef HAVE_MOAB_MESH
  else if (mesh_class == "MOAB")  
    {
      
      Teuchos::ParameterList moab_mesh_parameter_list = 
      	mesh_parameter_list.sublist("MOAB Mesh Parameters");
      
      string filename = moab_mesh_parameter_list.get<string>("Exodus file name");

      // shut up moab library..
      std::streambuf *store_buf = std::cout.rdbuf();
      std::cout.rdbuf(0);

      mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MOAB(filename.c_str(), MPI_COMM_WORLD));            
      std::cout.rdbuf(store_buf);
    }
#endif
#ifdef HAVE_STK_MESH
  else if (mesh_class == "STK")
    {
      string filename = mesh_parameter_list.get<string>("STK File name");
      
      mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_STK(MPI_COMM_WORLD, filename.c_str()));            
    }
#endif
  else
    {
      ostringstream error_message;
      error_message << "AMANZI_DEMO_DRIVER: main(): "
                    << "could not find mesh class \'" << mesh_class
                    << "\' is not enabled." << std::endl;
      Exceptions::amanzi_throw(Errors::Message(error_message.str()));
    }

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  // create the MPC
  Amanzi::MPC mpc(driver_parameter_list, mesh, obs_data);
  
  mpc.cycle_driver();
  
  delete comm;
      
}


