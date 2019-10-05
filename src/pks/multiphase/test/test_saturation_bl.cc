/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Quan Bui (mquanbui@math.umd.edu)
  Test for solving Buckley-Leverett equation with flux f(S) = S^2/(S^2 + (1-S)^2)
  Domain x = [0,1]; Initial condition S = 0 for x in [0,1]
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Saturation_PK.hh"
#include "OperatorDefs.hh"


/* **************************************************************** */
TEST(MULTIPHASE_SATURATION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: steady state 2-phase immiscible flow, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/test_saturation_bl.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  ParameterList slist = plist.sublist("State");
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "Cycle Driver");
  std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector();

  // time integration list
  ParameterList ti_list = plist.sublist("Cycle Driver").sublist("time intervals").sublist("TI 0");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 0.1, 0.1, mesh_size[0], 1, 1, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(slist));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  /* modify the default state for the problem at hand */
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";
  std::vector<int> ndofs(2, 1);
  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;

  std::string passwd("state");
  // Create darcy flux for upwind
  if (!S->HasField("darcy_flux")) {
    S->RequireField("darcy_flux", passwd)->SetMesh(mesh)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  // create solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  // Create saturation PK
  Multiphase::Saturation_PK* SPK = new Multiphase::Saturation_PK(*pk_tree_list, parameter_list, S, soln);
  SPK->Setup();
  S->Setup();
  S->InitializeFields();

  // modify flux
  Epetra_MultiVector& darcy_flux = *S->GetFieldData("darcy_flux", passwd)->ViewComponent("face", false);
  Point vel(1.0, 0.0, 0.0);
  for (int f = 0; f < darcy_flux.MyLength(); f++) {
    darcy_flux[0][f] = vel * mesh->face_normal(f);
  }
  S->GetField("darcy_flux", passwd)->set_initialized();

  // Init the PK and solve for saturation
  double dT, dT_actual, t_sim, t_final;
  //dT = 0.8/mesh_size[0];
  dT_actual = 0.0;
  //t_final = 0.6;
  dT = time_list->get<double>("time step");
  t_final = time_list->get<double>("end time");
  t_sim = 0.0;
  SPK->InitializeFields();
  S->CheckAllFieldsInitialized();
  Epetra_MultiVector u_0 = *S->GetFieldData("water_saturation", passwd)->ViewComponent("cell", false);
  SPK->InitializeSaturation();
  SPK->InitTimeInterval(ti_list);

  bool failed = true;
  while (!(t_sim > t_final || abs(t_sim - t_final) < 1e-8)) {
    failed = SPK->AdvanceStep(t_sim, t_sim + dT, false);
    if (!failed) {
      t_sim += dT;
      SPK->CommitState(S.ptr());
    } else {
      dT = dT / 2.0;
    }
  }
  std::cout << "Final time = " << t_sim << "\n";
  SPK->CommitState(S.ptr());

  // get the data for viz
  Epetra_MultiVector& ws_c = *S->GetFieldData("water_saturation", passwd)->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_saturation_bl.gmv");
    GMV::start_data();
    GMV::write_cell_data(u_0, 0, "initial_water_saturation");
    GMV::write_cell_data(ws_c, 0, "water_saturation");
    GMV::close_data_file();
  }

  FILE * f;
  f = fopen("ws.txt","w");
  for (int i = 0; i < ws_c.MyLength(); i++) {
    fprintf(f, "%e\n", ws_c[0][i]);
  }
  fflush(f);
  fclose(f);

  delete SPK;
}
