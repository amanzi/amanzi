/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

#include "NavierStokes_PK.hh"

/* **************************************************************** */
TEST(NAVIER_STOKES_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::NavierStokes;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D Navier Stokes" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/navier_stokes_2D.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  int nx = plist->get<int>("mesh resolution", 20);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, nx);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<NavierStokes_PK> NSPK = Teuchos::rcp(new NavierStokes_PK(plist, "navier stokes", S, soln));

  NSPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the Navier Stokes process kernel
  NSPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();
 
  // solve the problem 
  int itrs(0);
  int max_itrs = plist->get<int>("max iterations", 50);
  double T1 = plist->get<double>("end time", 100.0);
  double dT = plist->get<double>("initial time step", 1.0);
  double T(0.0), T0(0.0), dT0(dT), dTnext;

  // T = T1;
  while (T < T1 && itrs < max_itrs) {
    if (itrs == 0) {
      Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
      udot->PutScalar(0.0);
      NSPK->bdf1_dae()->SetInitialState(T0, soln, udot);

      NSPK->UpdatePreconditioner(T0, soln, dT0);
    }

    while (NSPK->bdf1_dae()->TimeStep(dT, dTnext, soln)) {
      dT = dTnext;
    }
    NSPK->bdf1_dae()->CommitSolution(dT, soln);

    T = NSPK->bdf1_dae()->time();
    dT = dTnext;
    itrs++;

    // reset primary fields
    Teuchos::RCP<PrimaryVariableFieldEvaluator> fluid_velocity_eval =
       Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator("pressure"));
    *S->GetFieldData("fluid_velocity", "navier stokes") = *soln->SubVector(0)->Data();
    fluid_velocity_eval->SetFieldAsChanged(S.ptr());

    Teuchos::RCP<PrimaryVariableFieldEvaluator> pressure_eval =
       Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator("pressure"));
    *S->GetFieldData("pressure", "navier stokes") = *soln->SubVector(1)->Data();
    pressure_eval->SetFieldAsChanged(S.ptr());

    // commit step
    NSPK->CommitStep(T - dT, T, S);
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  io.InitializeCycle(T, 1);
  const auto& u = *S->GetFieldData("fluid_velocity")->ViewComponent("node");
  const auto& p = *S->GetFieldData("pressure")->ViewComponent("cell");
  io.WriteVector(*u(0), "velocity_x", AmanziMesh::NODE);
  io.WriteVector(*u(1), "velocity_y", AmanziMesh::NODE);
  io.WriteVector(*p(0), "pressure", AmanziMesh::CELL);
  io.FinalizeCycle();
}
