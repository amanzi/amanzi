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
#include "IO.hh"
#include "GMVMesh.hh"
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* **************************************************************** */
void
RunTest(const std::string& filename, double tol)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D Richards, 2-layer model" << std::endl;

  // read parameter list
  std::string xmlFileName = filename;
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 6.0, 120.0, 3, 60);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree("flow");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(pk_tree, plist, S, soln));

  RPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // create the initial pressure function
  std::string passwd("");
  auto& pf = *S->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");
  auto& pm = *S->GetW<CompositeVector>("pressure_msp", passwd).ViewComponent("cell");

  for (int c = 0; c < pf.MyLength(); c++) {
    pm[0][c] = pf[0][c];
  }

  // initialize the Richards process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the problem
  TI_Specs ti_specs;
  ti_specs.T0 = 0.0;
  ti_specs.dT0 = 2000.0;
  ti_specs.T1 = 6.0e+10;
  ti_specs.max_itrs = 1000;

  AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  RPK->CommitStep(0.0, 1.0, Tags::DEFAULT); // dummy times

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string) "flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(pf, 0, "pressure");
    GMV::write_cell_data(pm, 0, "pressure_msp");
    GMV::close_data_file();
  }

  // check the pressure
  for (int c = 0; c < 60; c++) CHECK_CLOSE(pm[0][c], pf[0][c], tol);

  WriteStateStatistics(*S);
}


TEST(FLOW_2D_MULTISCALE_DPM)
{
  RunTest("test/flow_multiscale.xml", 0.2);
}

TEST(FLOW_2D_MULTISCALE_GDPM)
{
  RunTest("test/flow_multiscale_gdpm.xml", 0.8);
}
