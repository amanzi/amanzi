#include <iostream>
#include "stdlib.h"
#include "math.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "energy_tcm_registration.hh"
#include "energy_iem_registration.hh"
#include "eos_registration.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_registration.hh"
#include "pks_flow_registration.hh"
#include "pks_transport_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"


TEST(MPC_DRIVER_FLOW_MATRIX_FRACTURE) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  // read the main parameter list
  std::string xmlInFileName = "test/mpc_driver_single_fracture.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, &comm));

  // create mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory factory(&comm);
  factory.preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  // Create a verbose object to pass to the mesh_factory and mesh
  Teuchos::ParameterList mesh_params = plist -> sublist("mesh");
  // Read and initialize the unstructured mesh parameters
  Teuchos::ParameterList unstr_mesh_params = mesh_params.sublist("unstructured");
  // Decide on which mesh framework to use
  bool expert_params_specified = unstr_mesh_params.isSublist("expert");

  std::string file(""), format("");
  
  if (unstr_mesh_params.isSublist("read mesh file")) {
    Teuchos::ParameterList read_params = unstr_mesh_params.sublist("read mesh file");
    
    if (read_params.isParameter("file")) {
      file = read_params.get<std::string>("file");
    } else {
      std::cerr << "Must specify File parameter for Read option under mesh" << std::endl;
      throw std::exception();
    }

    if (read_params.isParameter("format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("format");

      if (format != "Exodus II" && format != "H5M") {	    
        std::cerr << "Can only read files in Exodus II or H5M format" << std::endl;
        throw std::exception();
      }
    } else {
      std::cerr << "Must specify 'format' parameter for Read option under mesh" << std::endl;
      throw std::exception();
    }

    int ierr;
    if (!file.empty()) {
      ierr = 0;
      try {
        // create the mesh from the file
        mesh = factory.create(file, gm);
	    
      } catch (const std::exception& e) {
        std::cerr << ": error: " << e.what() << std::endl;
        ierr++;
      }

    }
  }


  // create dummy observation data object
  Amanzi::ObservationData obs_data;    
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);
  
  // create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");

  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> mstk =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(mesh);
  Teuchos::RCP<AmanziMesh::Mesh> mesh_fracture =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(*mstk, names, AmanziMesh::FACE));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, &comm, obs_data);
  cycle_driver.Go();
}


