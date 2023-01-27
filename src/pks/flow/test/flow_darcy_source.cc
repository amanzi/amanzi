/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Tests with various sources: from a file
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
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"


/* *********************************************************************
* Turn on different sources and verify solution
********************************************************************* */
void
RunTestDarcySource(const std::string& xmlFileName)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D source model: " << xmlFileName << std::endl;

  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(-10.0, -5.0, 10.0, 0.0, 40, 10);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  auto soln = Teuchos::rcp(new TreeVector());
  auto DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("");
  auto& K = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");
  double diff_in_perm = 0.0;

  if (!S->GetRecord("permeability").initialized()) {
    for (int c = 0; c < K.MyLength(); c++) K[0][c] = 1.0;
    S->GetRecordW("permeability", "permeability").set_initialized();
  }

  // -- fluid density and viscosity
  S->GetW<double>("const_fluid_density", "state") = 1.0;
  S->GetRecordW("const_fluid_density", "state").set_initialized();

  S->GetW<double>("const_fluid_viscosity", "state") = 1.0;
  S->GetRecordW("const_fluid_viscosity", "state").set_initialized();

  // -- storativity
  S->GetW<CompositeVector>("specific_storage", passwd).PutScalar(0.1);
  S->GetRecordW("specific_storage", passwd).set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize();
  WriteStateStatistics(*S);

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 10; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    dt = DPK->get_dt();

    // verification
    double vmin, vmax;
    auto rhs =
      *DPK->my_operator(Operators::OperatorType::OPERATOR_MATRIX)->rhs()->ViewComponent("cell");
    rhs.MinValue(&vmin);
    rhs.MaxValue(&vmax);
    CHECK_CLOSE(vmin, vmax, 1e-12);
    if (n == 0) CHECK_CLOSE(vmin, -0.5, 1e-12);
  }
}


TEST(FLOW_2D_DARCY_SOURCE)
{
  RunTestDarcySource("test/flow_darcy_source.xml");
}
