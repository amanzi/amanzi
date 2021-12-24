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
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"

/* *********************************************************************
* Test of Darcy flow on polygonal mesh
********************************************************************* */
TEST(FLOW_2D_TRANSIENT_DARCY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: 2D transient Darcy, polygonal mesh" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_dual2D.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  pref.push_back(Framework::STK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create("test/dual2D.exo");

  // create a state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));

  DPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand 
  // -- permeability
  std::string passwd("flow"); 
  auto& K = *S->GetW<CompositeVector>("permeability", passwd).ViewComponent("cell");
  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 0.1;
    K[1][c] = 2.0;
  }
  S->GetRecordW("permeability", "flow").set_initialized();

  // -- fluid density and viscosity
  S->GetW<double>("const_fluid_density", passwd) = 1.0;
  S->GetRecordW("const_fluid_density", "flow").set_initialized();

  S->GetW<double>("const_fluid_viscosity", passwd) = 1.0;
  S->GetRecordW("const_fluid_viscosity", "flow").set_initialized();

  // -- gravity
  auto& gravity = S->GetW<AmanziGeometry::Point>("gravity", "state");
  gravity[1] = -1.0;
  S->GetRecordW("gravity", "state").set_initialized();

  // create the initial pressure function
  auto& p = *S->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = xc[1] * (xc[1] + 2.0);
  }
  S->GetRecordW("pressure", "flow").set_initialized();

  // initialize Darcy process kernel.
  DPK->Initialize(S.ptr());

auto vo = Teuchos::rcp(new Amanzi::VerboseObject("State", state_list));
WriteStateStatistics(*S, *vo);

  // transient solution
  double t_old(0.0), t_new, dt(0.1);
  for (int n = 0; n < 2; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, S);

    t_old = t_new;

    if (MyPID == 0) {
      GMV::open_data_file(*mesh, (std::string)"flow.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }
  }
}
