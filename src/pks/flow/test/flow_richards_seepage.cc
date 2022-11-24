/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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
#include "OutputXDMF.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

void
Flow2D_SeepageTest(std::string filename, bool deform)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D Richards, seepage boundary condition" << std::endl;

  // read parameter list and select left head
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(filename);

  // create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK, Framework::STK }));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 100.0, 50.0, 50, 25);

  // create an optional slop
  if (deform) {
    AmanziGeometry::Point xv(2);
    AmanziMesh::Entity_ID_List nodeids;
    AmanziGeometry::Point_List new_positions, final_positions;

    int nnodes = mesh->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
    for (int v = 0; v < nnodes; ++v) {
      mesh->node_get_coordinates(v, &xv);
      nodeids.push_back(v);
      xv[1] *= (xv[0] * 0.4 + (100.0 - xv[0])) / 100.0;
      new_positions.push_back(xv);
    }
    mesh->deform(nodeids, new_positions, false, &final_positions);
  }

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
  std::string passwd("");
  double rho = S->Get<double>("const_fluid_density");
  double g = (S->Get<AmanziGeometry::Point>("gravity"))[1];

  // create the initial pressure function
  auto& p = *S->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");

  double patm(101325.0), z0(30.0);
  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = patm + rho * g * (xc[1] - z0);
  }

  auto& lambda = *S->GetW<CompositeVector>("pressure", passwd).ViewComponent("face");
  RPK->DeriveFaceValuesFromCellValues(p, lambda);

  // create Richards process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the steady-state problem
  TI_Specs ti_specs;
  ti_specs.T0 = 0.0;
  ti_specs.dT0 = 10.0;
  ti_specs.T1 = 3000.0;
  ti_specs.max_itrs = 3000;

  AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  RPK->CommitStep(0.0, 1.0, Tags::DEFAULT); // dummy times for flow
  printf("seepage face total = %12.4f\n", RPK->seepage_mass());

  // output
  const auto& ws = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");

  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);
  io.InitializeCycle(ti_specs.T1, 1, "");
  io.WriteVector(*p(0), "pressure", AmanziMesh::CELL);
  io.WriteVector(*ws(0), "saturation", AmanziMesh::CELL);
  io.FinalizeCycle();

  /*
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::write_cell_data(ws, 0, "saturation");
    GMV::close_data_file();
  }
  */

  delete RPK;
}

TEST(FLOW_2D_RICHARDS_SEEPAGE)
{
  // Flow2D_SeepageTest("test/flow_richards_seepage_vertical.xml", false);
  Flow2D_SeepageTest("test/flow_richards_seepage.xml", true);
}
