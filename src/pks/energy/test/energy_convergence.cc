/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

// Flow
#include "EnergyTwoPhase_PK.hh"
#include "Analytic01.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Energy;


TEST(ENERGY_CONVERGENCE) {
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout <<"Convergence analysis on three random meshes" << std::endl;

  std::string xmlFileName = "test/energy_convergence.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // convergence estimate: Use n=3 for the full test.
  int nmeshes = plist->get<int>("number of meshes", 1);
  std::vector<double> h, p_error, v_error;

  for (int n = 0; n < nmeshes; n++) {
    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(2, region_list, comm);
    
    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);
    pref.push_back(STKMESH);

    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
    Teuchos::RCP<const Mesh> mesh;
    if (n == 0) {
      mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);
      // mesh = meshfactory("test/random_mesh1.exo", gm);
    } else if (n == 1) {
      mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 20, 20, gm);
      // mesh = meshfactory("test/random_mesh2.exo", gm);
    } else if (n == 2) {
      mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 40, 40, gm);
      // mesh = meshfactory("test/random_mesh3.exo", gm);
    }

    // create a simple state and populate it
    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("State");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::ParameterList pk_tree;
    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Teuchos::RCP<EnergyTwoPhase_PK> EPK = Teuchos::rcp(new EnergyTwoPhase_PK(pk_tree, plist, S, soln));

    EPK->Setup();
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    EPK->Initialize();
    S->CheckAllFieldsInitialized();

    // constant time stepping 
    int itrs(0);
    double t(0.0), t1(1.0), dt(0.1), dt_next;
    while (t < t1) {
      if (itrs == 0) {
        Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
        udot->PutScalar(0.0);
        EPK->bdf1_dae()->SetInitialState(t, soln, udot);
        EPK->UpdatePreconditioner(t, soln, dt);
      }

      EPK->bdf1_dae()->TimeStep(dt, dt_next, soln);
      CHECK(dt_next >= dt);
      EPK->bdf1_dae()->CommitSolution(dt, soln);

      t += dt;
      itrs++;
    }

    EPK->CommitStep(0.0, 1.0);

    // calculate errors
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
    Analytic01 ana(temp, mesh);

    double l2_norm, l2_err, inf_err;  // error checks
    ana.ComputeCellError(*temp->ViewComponent("cell"), t1, l2_norm, l2_err, inf_err);

    printf("mesh=%d bdf1_steps=%d  L2_temperature_err=%7.3e\n", n, itrs, l2_err);
    CHECK(l2_err < 2e-1);

    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(*temp->ViewComponent("cell"), 0, "temperature");
    GMV::close_data_file();
  }

  delete comm;
}

