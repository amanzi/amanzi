/*
  MultiPhase

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
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
#include "MeshFactory.hh"
#include "Mesh.hh"

#include "State.hh"
#include "OperatorDefs.hh"

// Multiphase
#include "MultiphaseReduced_PK.hh"


/* **************************************************************** */
TEST(MULTIPHASE_MODEL_I) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Multiphase;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase pk, model I" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/multiphase_model1.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 50, 5);

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Mumtiphase_PK", *plist));

  // create a simple state populate it
  auto state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  // create a solution vector
  ParameterList pk_tree = plist->sublist("PKs").sublist("multiphase");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  auto MPK = Teuchos::rcp(new MultiphaseReduced_PK(pk_tree, plist, S, soln));

  MPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the multiphase process kernel
  MPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();
  S->WriteStatistics(vo);

  // loop
  bool failed = true;
  double t(0.0), tend(1.0e+10), dt(tend / 10);
  while (t - tend) {
    bool failed = MPK->AdvanceStep(t, t + dt, false);

    t += dt;
    MPK->CommitStep(t, t + dt, S); 
    S->advance_time(dt);
    S->advance_cycle();

    if (MyPID == 0) {
      std::cout << "State time=" << S->time() << ", cycle=" << S->cycle() << std::endl;
    }
  }

  S->WriteStatistics(vo);
}
