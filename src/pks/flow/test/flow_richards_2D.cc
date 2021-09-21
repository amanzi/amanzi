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
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* **************************************************************** */
TEST(FLOW_2D_RICHARDS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D Richards, 2-layer model" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_richards_2D.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, -2.0, 1.0, 0.0, 18, 18);

  int itrs[2];
  for (int loop = 0; loop < 2; ++loop) {
    // create a simple state and populate it
    Teuchos::ParameterList state_list = plist->sublist("state");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(plist, "flow", S, soln));

    RPK->Setup(S.ptr());
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    // modify the default state for the problem at hand
    std::string passwd("flow"); 
    Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell");
  
    AmanziMesh::Entity_ID_List block;
    mesh->get_set_entities("Material 1", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);
    for (int i = 0; i != block.size(); ++i) {
      int c = block[i];
      K[0][c] = 0.1;
      K[1][c] = 2.0;
    }

    mesh->get_set_entities("Material 2", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);
    for (int i = 0; i != block.size(); ++i) {
      int c = block[i];
      K[0][c] = 0.5;
      K[1][c] = 0.5;
    }
    S->GetField("permeability", "flow")->set_initialized();

    // -- fluid vicosity
    S->GetFieldData("viscosity_liquid", passwd)->PutScalar(1.0);
    S->GetField("viscosity_liquid", passwd)->set_initialized();

    Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", "state");
    gravity[1] = -9.8;
    S->GetField("gravity", "state")->set_initialized();

    // create the initial pressure function
    Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell");

    for (int c = 0; c < p.MyLength(); c++) {
      const Point& xc = mesh->cell_centroid(c);
      p[0][c] = xc[1] * (xc[1] + 2.0);
    }

    // initialize the Richards process kernel
    RPK->Initialize(S.ptr());
    S->CheckAllFieldsInitialized();

    // solve the problem 
    TI_Specs ti_specs;
    ti_specs.T0 = 0.0;
    ti_specs.dT0 = 1.0;
    ti_specs.T1 = 100.0;
    ti_specs.max_itrs = 400;

    AdvanceToSteadyState(S, *RPK, ti_specs, soln);
    RPK->CommitStep(0.0, 1.0, S);  // dummy times
    itrs[loop] = ti_specs.num_itrs;

    if (MyPID == 0 && loop == 0) {
      GMV::open_data_file(*mesh, (std::string)"flow.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }

    // check the pressure
    int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c < ncells; c++) CHECK(p[0][c] > -4.0 && p[0][c] < 0.01);

    // modify the preconditioner
    plist->sublist("PKs").sublist("flow")
          .sublist("operators").sublist("diffusion operator").sublist("preconditioner")
          .set<std::string>("Newton correction", "approximate Jacobian");
  }

  // verify positive impact of Newton correction in the preconditioner.
  CHECK(itrs[1] < 0.7 * itrs[0]);
}
