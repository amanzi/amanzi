#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"

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

#include "MeshFactory.hh"
#include "State.hpp"
#include "MPC.hpp"

#include "errors.hh"
#include "exceptions.hh"

#include "amanzi_unstructured_grid_simulation_driver.hpp"
#include "InputParser.H"
#include "Teuchos_StrUtils.hpp"

Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run (const MPI_Comm&               mpi_comm,
                                             const Teuchos::ParameterList& input_parameter_list,
                                             Amanzi::ObservationData&      output_observations)
{
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // make sure only PE0 can write to std::cout
  int rank, ierr, aerr;
  MPI_Comm_rank(mpi_comm,&rank);

  if (rank!=0) {
    cout.rdbuf(0);
  } 

  ParameterList params_copy = Amanzi::AmanziInput::translate_state_sublist(input_parameter_list);

  // print parameter list
  std::cout << "======================> dumping parameter list <======================" << std::endl;
  Teuchos::writeParameterListToXmlOStream(params_copy, std::cout);
  std::cout << "======================> done dumping parameter list. <================"<<std::endl;

  using namespace std;

  Amanzi::AmanziMesh::MeshFactory factory(*comm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // get the Mesh sublist

  ierr = 0;
  Teuchos::ParameterList mesh_parameter_list = params_copy.sublist("Mesh");

  try {
    std::string framework = mesh_parameter_list.get<string>("Framework");
    Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::Simple)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::Simple);
    } else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MOAB)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MOAB);
    } else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::STKMESH)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::STKMESH);
    } else if (framework == Amanzi::AmanziMesh::framework_name(Amanzi::AmanziMesh::MSTK)) {
      prefs.clear(); prefs.push_back(Amanzi::AmanziMesh::MSTK);
    } else if (framework == "") {
      // do nothing
    } else {
      std::string s(framework);
      s += ": specified mesh framework preference not understood";
      amanzi_throw(Errors::Message(s));
    }
    factory.preference(prefs);
  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
    // do nothing, this means that the "Framework" parameter was not in the input
  } catch (const std::exception& e) {
    std::cerr << rank << ": error: " << e.what() << std::endl;
    ierr++;
  }

  comm->SumAll(&ierr, &aerr, 1);
  if (aerr > 0) {
    return Amanzi::Simulator::FAIL;
  }


  std::string file("");
  try {
    file = mesh_parameter_list.get<string>("Read");
  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
    // do nothing, this means that the "Read" parameter was not there
  }

  if (!file.empty()) {

                                // make a mesh from a mesh file

    ierr = 0;
    try {
      mesh = factory.create(file);
    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }
  
    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) {
      return Amanzi::Simulator::FAIL;
    }
  
  } else {

                                // generate a hex mesh

    ierr = 0;

    try {
      mesh = factory(mesh_parameter_list.sublist("Generate"));
    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }
  
    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) {
      return Amanzi::Simulator::FAIL;
    }

  }

  ASSERT(!mesh.is_null());

  // create the MPC
  Amanzi::MPC mpc(params_copy, mesh, output_observations);
  
  mpc.cycle_driver();
  
  mesh.reset();
  delete comm;
      
  return Amanzi::Simulator::SUCCESS;
}


