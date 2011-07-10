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

#include "amanzi_unstructured_grid_simulation_driver.hpp"

Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run (const MPI_Comm&               mpi_comm,
                                             const Teuchos::ParameterList& input_parameter_list,
                                             ObservationData&              output_observations)
{
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // make sure only PE0 can write to std::cout
  int rank;
  MPI_Comm_rank(mpi_comm,&rank);

  if (rank!=0) {
    cout.rdbuf(0);
  } 

  // print parameter list
  std::cout << "======================> dumping parameter list <======================" << std::endl;
  Teuchos::writeParameterListToXmlOStream(input_parameter_list, std::cout);
  std::cout << "======================> done dumping parameter list. <================"<<std::endl;

  using namespace std;

  // get the Mesh sublist
  Teuchos::ParameterList mesh_parameter_list = input_parameter_list.sublist("Mesh");

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


  // create the MPC
  Amanzi::MPC mpc(input_parameter_list, mesh);
  
  mpc.cycle_driver();
  
  delete comm;
      
}


