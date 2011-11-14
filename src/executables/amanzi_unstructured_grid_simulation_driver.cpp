#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "State.hpp"
#include "MPC.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "amanzi_unstructured_grid_simulation_driver.hpp"
#include "InputParserIS.H"

Amanzi::Simulator::ReturnType
AmanziUnstructuredGridSimulationDriver::Run (const MPI_Comm&               mpi_comm,
                                             Teuchos::ParameterList& input_parameter_list,
                                             Amanzi::ObservationData&      output_observations)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab


#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(mpi_comm);
#else  
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int rank, ierr, aerr;
  MPI_Comm_rank(mpi_comm,&rank);

  bool native = input_parameter_list.get<bool>("Native Unstructured Input",false);
  
  Teuchos::ParameterList new_list; 
  Teuchos::ParameterList sub_list;
  
  if (! native)
    {

      new_list = Amanzi::AmanziInput::translate(input_parameter_list);

      // Amanzi::AmanziInput::init_global_info(input_parameter_list);

      // sub_list = Amanzi::AmanziInput::create_Checkpoint_Data_List (input_parameter_list);
      // new_list.sublist("Checkpoint Data") = sub_list;    
      
      // sub_list = Amanzi::AmanziInput::create_Visualization_Data_List (input_parameter_list);
      // new_list.sublist("Visualization Data") = sub_list;    

      // sub_list = Amanzi::AmanziInput::create_Observation_Data_List (input_parameter_list); 
      // new_list.sublist("Observation Data") = sub_list;

      // sub_list = Amanzi::AmanziInput::get_Regions_List(input_parameter_list);
      // new_list.sublist("Regions") = sub_list;
      
      // sub_list = Amanzi::AmanziInput::get_Mesh_List(input_parameter_list);
      // new_list.sublist("Mesh") = sub_list;

      // sub_list = Amanzi::AmanziInput::get_Domain_List(input_parameter_list);
      // new_list.sublist("Domain") = sub_list;

      // sub_list = Amanzi::AmanziInput::create_MPC_List(input_parameter_list);
      // new_list.sublist("MPC") = sub_list;
      
      // sub_list = Amanzi::AmanziInput::create_Transport_List(input_parameter_list);
      // new_list.sublist("Transport") = sub_list;
       
      // sub_list = Amanzi::AmanziInput::create_State_List(input_parameter_list);
      // new_list.sublist("State") = sub_list;

      // sub_list = Amanzi::AmanziInput::create_Flow_List(input_parameter_list);
      // new_list.sublist("Flow") = sub_list;

      //params_copy = Amanzi::AmanziInput::translate_state_sublist(input_parameter_list);
    }
  else
    {
      new_list = input_parameter_list;
    }

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))	  
    {  
      // print parameter list
      *out << "======================> dumping parameter list <======================" << std::endl;
      Teuchos::writeParameterListToXmlOStream(new_list, *out);
      *out << "======================> done dumping parameter list. <================"<<std::endl;
    }
  return Amanzi::Simulator::SUCCESS;
  using namespace std;

  Amanzi::AmanziMesh::MeshFactory factory(*comm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // get the Mesh sublist

  ierr = 0;
  Teuchos::ParameterList mesh_parameter_list = new_list.sublist("Mesh");

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
  Amanzi::MPC mpc(new_list, mesh, comm, output_observations);
  
  mpc.cycle_driver();
  
  mesh.reset();
  delete comm;
      
  return Amanzi::Simulator::SUCCESS;
}


