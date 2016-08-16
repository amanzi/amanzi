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
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"


/* *********************************************************************
* Two tests with different time step controllers.
********************************************************************* */
void RunTestMarshak(std::string controller) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "\nTest: 2D well model: " << controller << std::endl;

  // read parameter list
  std::string xmlFileName = controller;
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(-10.0, -5.0, 10.0, 0.0, 101, 50, gm);

  // create a simple state and populate it
  Amanzi::VerboseObject::hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S));
  DPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("flow"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  
  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 0.1;
    K[1][c] = 2.0;
  }
  S->GetField("permeability", "flow")->set_initialized();

  // -- fluid density and viscosity
  *S->GetScalarData("fluid_density", passwd) = 1.0;
  S->GetField("fluid_density", "flow")->set_initialized();

  *S->GetScalarData("fluid_viscosity", passwd) = 1.0;
  S->GetField("fluid_viscosity", "flow")->set_initialized();

  // -- gravity
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", "state");
  gravity[1] = -1.0;
  S->GetField("gravity", "state")->set_initialized();

  // -- storativity
  S->GetFieldData("specific_storage", passwd)->PutScalar(0.1);
  S->GetField("specific_storage", "flow")->set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize(S.ptr());

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 10; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, S);

    t_old = t_new;

    if (MyPID == 0) {
      const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell");
      GMV::open_data_file(*mesh, (std::string)"flow2D.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }

    dt = DPK->get_dt();
std::cout << "dT=" << dt << std::endl;
  }
}


TEST(FLOW_2D_DARCY_WELL_STANDARD) {
  RunTestMarshak("test/flow_darcy_well.xml");
}

TEST(FLOW_2D_DARCY_WELL_ADAPRIVE) {
  RunTestMarshak("test/flow_darcy_well_adaptive.xml");
}


/* **************************************************************** */
TEST(FLOW_3D_DARCY_WELL) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "Test: 3D Darcy flow, two wells" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_well_3D.xml";
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an MSTK mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(-10.0, -1.0, -5.0, 10.0, 1.0, 0.0, 101, 10, 25, gm);

  // create a simple state and populate it
  Amanzi::VerboseObject::hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S));
  DPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("flow"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  
  for (int c = 0; c < K.MyLength(); c++) {
    K[0][c] = 0.1;
    K[1][c] = 2.0;
    K[2][c] = 2.0;
  }
  S->GetField("permeability", "flow")->set_initialized();

  // -- fluid density and viscosity
  *S->GetScalarData("fluid_density", passwd) = 1.0;
  S->GetField("fluid_density", "flow")->set_initialized();

  *S->GetScalarData("fluid_viscosity", passwd) = 1.0;
  S->GetField("fluid_viscosity", "flow")->set_initialized();

  // -- gravity
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", "state");
  gravity[2] = -1.0;
  S->GetField("gravity", "state")->set_initialized();

  // -- storativity
  S->GetFieldData("specific_storage", passwd)->PutScalar(0.1);
  S->GetField("specific_storage", "flow")->set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize(S.ptr());

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 2; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, S);

    t_old = t_new;

    if (MyPID == 0) {
      const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell");
      GMV::open_data_file(*mesh, (std::string)"flow3D.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }
  }
}
