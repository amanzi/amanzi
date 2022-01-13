/*
  Energy

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
#include "UnitTest++.h"

// Amanzi
#include "CompositeVector.hh"
#include "IO.hh"
#include "EnergyOnePhase_PK.hh"
#include "MeshFactory.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "State.hh"
#include "VerboseObject.hh"
#include "WhetStoneDefs.hh"

/* **************************************************************** 
* Runs to a steady state
* ************************************************************** */
TEST(ENERGY_ONE_PHASE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;
  using namespace Amanzi::Energy;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: steady state calculation" << std::endl;

  // read parameter list 
  std::string xmlFileName = "test/energy_one_phase.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::rcp(new Teuchos::ParameterList(xmlreader.getParameters()));

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 10, 2);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("energy");
  auto soln = Teuchos::rcp(new TreeVector());
  auto EPK = Teuchos::rcp(new EnergyOnePhase_PK(pk_tree, plist, S, soln));

  EPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  EPK->Initialize();
  S->CheckAllFieldsInitialized();

  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("EnergyOnePhase", *plist));
  WriteStateStatistics(*S, *vo);

  // constant time stepping 
  int itrs(0);
  double t(0.0), dt(0.1), t1(5.5), dt_next;
  while (t < t1) {
    // swap conserved quntity (no backup, we check dt_next instead)
    const auto& e = S->Get<CompositeVector>("energy");
    auto& e_prev = S->GetW<CompositeVector>("prev_energy", "thermal");
    e_prev = e;

    if (itrs == 0) {
      auto udot = Teuchos::rcp(new TreeVector(*soln));
      udot->PutScalar(0.0);
      EPK->bdf1_dae()->SetInitialState(t, soln, udot);
      EPK->UpdatePreconditioner(t, soln, dt);
    }

    EPK->bdf1_dae()->TimeStep(dt, dt_next, soln);
    CHECK(dt_next >= dt);
    EPK->bdf1_dae()->CommitSolution(dt, soln);
    Teuchos::rcp_static_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace> >(
        S->GetEvaluatorPtr("temperature"))->SetChanged();

    t += dt;
    itrs++;
  }

  EPK->CommitStep(0.0, 1.0, Tags::DEFAULT);
  WriteStateStatistics(*S, *vo);

  auto temp = *S->Get<CompositeVector>("temperature").ViewComponent("cell");
  for (int c = 0; c < 20; ++c) { 
    CHECK_CLOSE(1.5, temp[0][c], 2e-8);
  }
}
