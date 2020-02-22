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
#include "OutputXDMF.hh"

// Multiphase
#include "MultiphaseModelI_PK.hh"


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
  std::string xmlFileName = "test/multiphase_modelIa.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a MSTK mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  // RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 200, 5);
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 50, 3);

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Multiphase_PK", *plist));

  // create a simple state populate it
  auto state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  // create a solution vector
  ParameterList pk_tree = plist->sublist("PKs").sublist("multiphase");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  auto MPK = Teuchos::rcp(new MultiphaseModelI_PK(pk_tree, plist, S, soln));

  MPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the multiphase process kernel
  MPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();
  S->WriteStatistics(vo);

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh, true, false));

  // loop
  int iloop(0);
  bool failed = true;
  double t(0.0), tend(1.57e+12), dt(1.5768e+7), dt_max(3e+10);
  while (t < tend && iloop < 400) {
    while (MPK->AdvanceStep(t, t + dt, false)) { dt /= 4; }

    MPK->CommitStep(t, t + dt, S);
    S->advance_cycle();

    S->advance_time(dt);
    t += dt;
    dt = std::min(dt_max, dt * 1.2);
    iloop++;

    // output solution
    if (iloop % 5 == 0) {
      io->InitializeCycle(t, iloop);
      const auto& u0 = *S->GetFieldData("pressure_liquid")->ViewComponent("cell");
      const auto& u1 = *S->GetFieldData("saturation_liquid")->ViewComponent("cell");
      const auto& u2 = *S->GetFieldData("mole_fraction_gas")->ViewComponent("cell");

      io->WriteVector(*u0(0), "pressure", AmanziMesh::CELL);
      io->WriteVector(*u1(0), "saturation", AmanziMesh::CELL);
      io->WriteVector(*u2(0), "mole fraction gas", AmanziMesh::CELL);
      io->FinalizeCycle();

      S->WriteStatistics(vo);
    }
  }

  S->WriteStatistics(vo);

  // verification
  double dmin, dmax;
  const auto& sl = *S->GetFieldData("saturation_liquid")->ViewComponent("cell");
  sl.MinValue(&dmin);
  sl.MaxValue(&dmax);
  CHECK(dmin >= 0.0 && dmax <= 1.0);
  
  S->GetFieldData("ncp_fg")->NormInf(&dmax);
  CHECK(dmax <= 1.0e-14);

  const auto& xg = *S->GetFieldData("mole_fraction_gas")->ViewComponent("cell");
  xg.MinValue(&dmin);
  xg.MaxValue(&dmax);
  CHECK(dmin >= 0.0 && dmax <= 1.0);
}
