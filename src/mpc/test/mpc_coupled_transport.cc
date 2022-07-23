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
#include "eos_registration.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_transport_registration.hh"
#include "State.hh"


TEST(MPC_DRIVER_TRANSPORT_MATRIX_FRACTURE) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  
  std::string xmlInFileName = "test/mpc_coupled_transport.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({Framework::MSTK}));
  auto mesh = factory.create(0.0, 0.0, 0.0, 6.0, 6.0, 6.0, 6, 6, 6);

  // create dummy observation data object
  Amanzi::ObservationData obs_data;    
  
  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  //create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture = factory.create(mesh, names, AmanziMesh::FACE);

  S->RegisterMesh("fracture", mesh_fracture);

  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();

  const auto& tcc_m = *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  double cmin(1e+99);
  double cmax = -1e+99;
  for (int c = 0; c < 216; ++c) {
    const auto& xc = mesh->cell_centroid(c);
    if (std::fabs(xc[0] - 5.5) < 1e-3) {
      cmin = std::min(cmin, tcc_m[0][c]);
      cmax = std::max(cmin, tcc_m[0][c]);
    }
  }
  CHECK_CLOSE(cmin, cmax, 1e-12);

  const auto& tcc_f = *S->Get<CompositeVector>("fracture-total_component_concentration").ViewComponent("cell");
  for (int c = 0; c < 36; ++c) {
    const auto& xc = mesh->cell_centroid(c);
    if (std::fabs(xc[0] - 5.5) < 1e-3) {
      cmin = std::min(cmin, tcc_f[0][c]);
      cmax = std::max(cmin, tcc_f[0][c]);
    }
  }
  CHECK_CLOSE(cmin, cmax, 1e-12);
}


