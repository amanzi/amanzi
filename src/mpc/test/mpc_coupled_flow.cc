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
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
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

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  
  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_coupled_flow.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({Framework::MSTK}));
  auto mesh = factory.create(0.0, 0.0, 0.0, 216.0, 10.0, 120.0, 9, 2, 20, true, true);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;    
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  // auto mesh_fracture = factory.create(mesh, names, AmanziMesh::FACE);
  auto mesh_fracture = Teuchos::rcp(new MeshExtractedManifold(
      mesh, "fracture", AmanziMesh::FACE, comm, gm, mesh_list, true, false));

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  // test pressure in fracture (5% error)
  double p0 = 101325.0;
  double rho = 998.2;
  double mu = 0.001002;
  double K1 = 1.0e-11;
  double kn = 4.0e-8;
  double L = 120.0;
  double q0 = -2.0e-3;
  double gravity = 9.81;

  // test pressure in fracture
  K1 *= rho / mu;
  double pf_exact = p0 - q0 * (L / K1 / 2 + 1.0 / kn) - gravity * rho * L / 2;

  auto& pf = *S->GetFieldData("fracture-pressure")->ViewComponent("cell");
  for (int c = 0; c < pf.MyLength(); ++c) {
    if (c == 0) std::cout << "Fracture pressure: " << pf[0][c] << ",  exact: " << pf_exact << std::endl;
    CHECK(std::fabs(pf[0][c] - pf_exact) < 1e-8 * std::fabs(pf_exact));
  }

  // test flux in bottom domain
  Amanzi::AmanziGeometry::Point uf_exact(0.0, 0.0, q0);

  auto& uf = *S->GetFieldData("darcy_flux")->ViewComponent("face");
  const auto& fmap = *S->GetFieldData("darcy_flux")->ComponentMap("face");

  bool flag(true);
  int nfaces = mesh->num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int f = 0; f < nfaces; ++f) {
    const auto& normal = mesh->face_normal(f);

    int g = fmap.FirstPointInElement(f);
    double flux = (uf_exact * normal) / rho;

    if (fmap.ElementSize(f) == 2 && flag) {
       std::cout << "Matrix Darcy fluxes: " << uf[0][g] << ",  exact: " << flux << std::endl;
       flag = false;
    }
    CHECK(std::fabs(uf[0][g] - flux) < 1e-8 * std::fabs(q0) + 1e-14);
  }
}


