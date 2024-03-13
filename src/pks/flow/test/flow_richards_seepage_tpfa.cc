/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* **************************************************************** */
TEST(FLOW_2D_RICHARDS_SEEPAGE_TPFA)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D Richards(FV), seepage boundary condition" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_richards_seepage_tpfa.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 100.0, 50.0, 100, 50);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Richards_PK* RPK = new Richards_PK(plist, "flow", S, soln);

  RPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("");
  auto& K = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");

  for (int c = 0; c != K.MyLength(); ++c) {
    K[0][c] = 5.0e-13;
    K[1][c] = 5.0e-14;
  }
  S->GetRecordW("permeability", "permeability").set_initialized();

  // -- fluid density and viscosity
  double rho = S->Get<double>("const_fluid_density");

  S->GetW<CompositeVector>("viscosity_liquid", "viscosity_liquid").PutScalar(0.00089);

  // -- gravity
  auto gravity = S->Get<AmanziGeometry::Point>("gravity");
  double g = gravity[1];

  // create the initial pressure function
  auto& p = *S->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");

  double p0(101325.0), z0(30.0);
  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->getCellCentroid(c);
    p[0][c] = p0 + rho * g * (xc[1] - z0);
  }

  // create Richards process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the steady-state problem
  TI_Specs ti_specs;
  ti_specs.T0 = 0.0;
  ti_specs.dT0 = 1.0;
  ti_specs.T1 = 1e+7;
  ti_specs.max_itrs = 10;

  AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  RPK->set_dt(1.0);
  printf("time step = %12.4f\n", RPK->get_dt());
  printf("seepage face total = %12.4f\n", RPK->seepage_mass());
  RPK->CommitStep(0.0, 1.0, Tags::DEFAULT); // dummy times for flow

  const auto& ws = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string) "flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::write_cell_data(ws, 0, "saturation");
    GMV::close_data_file();
  }

  delete RPK;
}
