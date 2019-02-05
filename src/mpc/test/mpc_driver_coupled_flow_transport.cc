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
  std::string xmlInFileName = "test/mpc_driver_coupled_flow_transport.xml";
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
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = factory(0.0, 0.0, 0.0, 216.0, 10.0, 120.0, 3, 2, 10, gm);
  // Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = factory("test/single_fracture_tet.exo", gm);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;    
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  Amanzi::MeshAudit mesh_auditor(mesh);
  int status = mesh_auditor.Verify();
  if (status == 0) {
    std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
  } else {
    Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
    Exceptions::amanzi_throw(msg);
  }
  
  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");

  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> mstk =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(mesh);
  Teuchos::RCP<AmanziMesh::Mesh> mesh_fracture =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(*mstk, names, AmanziMesh::FACE));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, &comm, obs_data);
  cycle_driver.Go();

  // test pressure in fracture (5% error)
  double p0 = 101325.0;
  double rho = 998.2;
  double mu = 0.001002;
  double K1 = 1.0e-11;
  double kn = 4.0e-8;
  double L = 60.0;
  double q0 = -2.0e-3;

  K1 *= rho / mu;
  double pf_exact = p0 - q0 * (L / K1 / 2 + 1.0 / kn);

  // double pf = (*S->GetFieldData("fracture-pressure")->ViewComponent("cell"))[0][0];
  // std::cout << "Fracture pressure: " << pf << ",  exact: " << pf_exact << std::endl;
  // CHECK(std::fabs(pf - pf_exact) < 0.05 * std::fabs(pf_exact));
}


