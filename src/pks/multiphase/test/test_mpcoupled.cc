/*
  This is the flow component of the Amanzi code. 

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

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "MPCoupled_PK.hh"
#include "OperatorDefs.hh"


/* **************************************************************** */
TEST(MULTIPHASE_MPCOUPLED) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase coupled, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  //std::string xmlFileName = "test/test_pressure_pk.xml";
  //std::string xmlFileName = "test/test_mpcoupled.xml";
  std::string xmlFileName = "test/test_grav_1d.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "Cycle Driver");
  ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  std::vector<double> high_coor = mesh_list.get<Teuchos::Array<double> > ("Domain High Coordinate").toVector();
  std::vector<double> low_coor = mesh_list.get<Teuchos::Array<double> > ("Domain Low Coordinate").toVector();
  std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector();

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(low_coor[0], low_coor[1],
                                     high_coor[0], high_coor[1],
                                     mesh_size[0], mesh_size[1], gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(*state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  std::string passwd("state");
  
  // create solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  // Create pressure PK
  Multiphase::MPCoupled_PK* MPK = new Multiphase::MPCoupled_PK(*pk_tree_list, parameter_list, S, soln);
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  MPK->Initialize();

  double dT, t_sim, t_final;
  bool failed = true;
  //double dT = 0.05;
  //double dT = 0.6/mesh_size[0];
  //dT = 0.1;
  t_sim = 0.0;
  //t_final = 0.2;
  dT = time_list->get<double>("time step");
  t_final = time_list->get<double>("end time");
  while (!(t_sim > t_final || abs(t_sim - t_final) < 1e-8)) {
    failed = MPK->AdvanceStep(t_sim, t_sim + dT, false);
    if (failed) {
      dT /= 2.0;
    } else {
      t_sim += dT;
      MPK->CommitStep(t_sim, t_sim + dT); 
    } 
  }
  std::cout << "t_sim: " << t_sim << "\n"; 

  // get the data for viz
  Epetra_MultiVector& p1 = *S->GetFieldData("phase1_pressure",passwd)->ViewComponent("cell");
  Epetra_MultiVector& p2 = *S->GetFieldData("phase2_pressure",passwd)->ViewComponent("cell");
  Epetra_MultiVector& s1 = *S->GetFieldData("phase1_saturation",passwd)->ViewComponent("cell");
  Epetra_MultiVector& s2 = *S->GetFieldData("phase2_saturation",passwd)->ViewComponent("cell");
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_mpcoupled.gmv");
    GMV::start_data();
    GMV::write_cell_data(p1, 0, "p1");
    GMV::write_cell_data(s1, 0, "s1");  
    GMV::write_cell_data(p2, 0, "p2");
    GMV::write_cell_data(s2, 0, "s2");    
    GMV::close_data_file();
  }

  FILE * out;
  out = fopen("ws.txt","w");
  for (int i = 0; i < s1.MyLength(); i++) {
    fprintf(out, "%e\n", s1[0][i]);
  }
  fflush(out);
  fclose(out);
  
  delete MPK;
}
