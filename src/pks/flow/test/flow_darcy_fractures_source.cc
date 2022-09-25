/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"

/* **************************************************************** */
TEST(DARCY_TWO_FRACTURES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Flow;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: Darcy PK in fractures: uniform water injection" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_fractures_source.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  // extract fractures mesh
  std::vector<std::string> setnames({"fracture 1" });
  RCP<Mesh> mesh = meshfactory.create(mesh3D, setnames, AmanziMesh::FACE);

  // create state and initialize
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("fracture", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup();

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  DPK->Initialize();
  S->CheckAllFieldsInitialized();
  WriteStateStatistics(*S);

  // time stepping
  double t_old(0.0), t_new, dt(10.0);
  for (int n = 0; n < 10; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, Tags::DEFAULT);
    S->set_time(t_new);

    t_old = t_new;
  }
  WriteStateStatistics(*S);

  // double V = 0.25;   // fracture area [m^2]
  // double a = 0.01;   // aprerture [m]
  double Q = 8.0e-2;    // source [kg/s]
  double Ss = 0.002; // specific storage [m^-1]
  double g = 10.0;   // gravity [m/s^2]
  double p_old = 200000.0;
  // double p_new = p_old + (dt * 10) * (Q / V) * g / (Ss * a); 
  double p_new = p_old + (dt * 10) * (Q) * g / (Ss); 

  auto& p = *S->GetW<CompositeVector>("fracture-pressure", Tags::DEFAULT, "flow").ViewComponent("cell");
  for (int c = 0; c < p.MyLength(); c++) {
    CHECK_CLOSE(p[0][c], p_new, 1.0);
  }
}





