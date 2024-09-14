/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cmath>

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "CycleDriver.hh"
#include "eos_reg.hh"
#include "evaluators_mpc_reg.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "models_flow_reg.hh"
#include "models_transport_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"


void
runTest()
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  std::cout << "\nTEST: copuled flow and transport, implicit scheme for gas" << std::endl;

  // read and modify input list
  std::string xmlInFileName = "test/mpc_coupled_flow_transport_gas.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, -5.0, 0.0, 100.0, 5.0, 100.0, 20, 2, 25);

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  // create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture = factory.create(mesh, names, AmanziMesh::Entity_kind::FACE);
  S->RegisterMesh("fracture", mesh_fracture);

  // create PK
  Teuchos::ParameterList tmp =
    plist->sublist("cycle driver").sublist("time periods").sublist("TP 0").sublist("PK tree");
  std::string name = tmp.begin()->first;
  Teuchos::ParameterList pk_tree = tmp.sublist(name);

  auto soln = Teuchos::rcp(new TreeVector());
  auto pk = Teuchos::rcp(new FlowReactiveTransport_PK(pk_tree, plist, S, soln));

  pk->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  pk->Initialize();
  S->CheckAllFieldsInitialized();
  pk->CommitStep(0.0, 1.0, Tags::DEFAULT);
  pk->CalculateDiagnostics(Tags::DEFAULT);

  // perform given steps
  double dt(1.0), told(0.0), tnew;
  double f1max, f2max;

  for (int n = 0; n < 100; ++n) {
    tnew = told + dt;
    pk->AdvanceStep(told, tnew, true);
    pk->CommitStep(told, tnew, Tags::DEFAULT);
    told = tnew;
    dt = std::min(7200.0, dt * 1.2);

    S->GetEvaluator("molal_concentration").Update(*S, "mpc");
    S->GetEvaluator("fracture-molal_concentration").Update(*S, "mpc");

    const auto& f1 = S->Get<CompositeVector>("molal_concentration", Tags::DEFAULT);
    const auto& f2 = S->Get<CompositeVector>("fracture-molal_concentration", Tags::DEFAULT);
    f1.NormInf(&f1max);
    f2.NormInf(&f2max);

    CHECK(f2max < f1max);
    std::cout << "Max molal conc: " << f1max << " " << f2max << "  t = " << tnew << std::endl;
  }

  WriteStateStatistics(*S);
}


TEST(MPC_COUPLED_FLOW_TRANSPORT_GAS)
{
  runTest();
}
