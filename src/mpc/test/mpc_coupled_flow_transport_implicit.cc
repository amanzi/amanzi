#include <iostream>
#include <cstdlib>
#include <cmath>

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "eos_registration.hh"
#include "mdm_transport_registration.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "mpc_pks_registration.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_registration.hh"
#include "pks_transport_registration.hh"
#include "State.hh"
#include "wrm_flow_registration.hh"


TEST(MPC_DRIVER_COUPLED_FLOW_TRANSPORT_IMPLICIT) {

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  
  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_coupled_flow_transport_implicit.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);
  
  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({Framework::MSTK}));
  auto mesh = factory.create(0.0, 0.0, 0.0, 216.0, 10.0, 120.0, 9, 2, 20);

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

  // test pressure in fracture (5% error)
  double p0 = 101325.0;
  double rho = 998.2;
  double mu = 0.001002;
  double K1 = 1.0e-11;
  double kn = 4.0e-8;
  double L = 120.0;
  double q0 = -2.0e-3;
  double g = 9.81;

  K1 *= rho / mu;
  double pf_exact = p0 - q0 * (L / K1 / 2 + 1.0 / kn) - g * rho * L / 2;

  double pf = (*S->Get<CompositeVector>("fracture-pressure").ViewComponent("cell"))[0][0];
  std::cout << "Fracture pressure: " << pf << ",  exact: " << pf_exact << std::endl;
  CHECK(std::fabs(pf - pf_exact) < 0.05 * std::fabs(pf_exact));

  // checking that we created only three PKs each time period
  CHECK(PKFactory::num_pks == 6);
  std::cout << PKFactory::list_pks << std::endl;
}


