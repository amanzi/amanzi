#include <iostream>
#include "stdlib.h"
#include "math.h"

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

#include "dbc.hh"
#include "errors.hh"
#include "MeshFactory.hh"
#include "State.hpp"
#include "MPC.hpp"
#include "Region.hh"
#include "RectangularRegion.hh"
#include "MeshDefs.hh"

#include "Teuchos_MPISession.hpp"
#include "InputParser.H"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);

  // make sure only PE0 can write to std::cout
  int rank, ierr, aerr;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // if (rank!=0) {
  //   cout.rdbuf(0);
  // }
 
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
  
  // read the main parameter list
  Teuchos::ParameterList driver_parameter_list;

  ierr = 0;
  try {
    Teuchos::updateParametersFromXmlFile(xmlInFileName,&driver_parameter_list);
  } catch (const std::runtime_error& e) {
    std::cerr << rank << ": error: " << e.what() << std::endl;
    ierr++;
  }

  comm->SumAll(&ierr, &aerr, 1);
  if (aerr > 0) {
    return 1;
  }
    
    
  ParameterList trans_parameter_list = Amanzi::AmanziInput::translate_state_sublist(driver_parameter_list);
  driver_parameter_list = trans_parameter_list;


  // print parameter list
  std::cout << "======================> dumping parameter list <======================" << std::endl;
  Teuchos::writeParameterListToXmlOStream(driver_parameter_list, std::cout);
  std::cout << "======================> done dumping parameter list. <================"<<std::endl;

  Amanzi::AmanziMesh::MeshFactory factory(comm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  
  // get the Mesh sublist

  ierr = 0;
  Teuchos::ParameterList mesh_parameter_list = driver_parameter_list.sublist("Mesh");

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
    return 3;
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
      return 3;
    }
  
  } else {

                                // generate a hex mesh

    ierr = 0;

    try {

      Teuchos::ParameterList& generate_parameter_list(mesh_parameter_list.sublist("Generate"));
      mesh = factory(generate_parameter_list);
    } catch (const std::exception& e) {
      std::cerr << rank << ": error: " << e.what() << std::endl;
      ierr++;
    }
  
    comm->SumAll(&ierr, &aerr, 1);
    if (aerr > 0) {
      return 3;
    }

  }

  ASSERT(!mesh.is_null());

#if 0
  std::map<std::string,Amanzi::AmanziMesh::Set_ID> region_label_ID_map;
  Amanzi::AmanziMesh::Entity_kind kind = Amanzi::AmanziMesh::CELL;
  for (RegionMap::const_iterator it=region_map.begin(); it!=region_map.end(); ++it) {

    // Does not exist yet
    // region_label_ID_map[ it->first ] = mesh->define_set( it->second, kind );
  }
#endif

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  // create the MPC
  Amanzi::MPC mpc(driver_parameter_list, mesh, comm, obs_data);
  
  mpc.cycle_driver();
  
  delete comm;
      
}


