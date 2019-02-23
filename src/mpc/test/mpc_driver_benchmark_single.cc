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
#include "MeshAudit.hh"
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
  
  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_driver_benchmark_single.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, &comm));

  // create mesh
  MeshFactory factory(&comm);

  factory.preference(FrameworkPreference({Framework::MSTK}));
  factory.set_partitioner(Amanzi::AmanziMesh::Partitioner_type::ZOLTAN_GRAPH);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = factory("test/single_fracture_ref2.exo", gm); 
  // std::string meshfile = plist->sublist("mesh").sublist("unstructured").sublist("read mesh file").get<std::string>("file");
  // Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = factory(meshfile, gm);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;    
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));

  S->RegisterMesh("domain", mesh);

  /*  
  Amanzi::MeshAudit mesh_auditor(mesh);
  int status = mesh_auditor.Verify();
  if (status != 0) {
    Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
    Exceptions::amanzi_throw(msg);
  }
  */
  
  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");

  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> mstk =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(mesh);
  Teuchos::RCP<AmanziMesh::Mesh> mesh_fracture =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(&*mstk, names, AmanziMesh::FACE));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, &comm, obs_data);
  cycle_driver.Go();
}


