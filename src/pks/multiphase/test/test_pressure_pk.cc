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
#include "Pressure_PK.hh"
#include "OperatorDefs.hh"

// Analytic data 2
// u(x,y) = 2x + 4y
double u_exact(const Amanzi::AmanziGeometry::Point& p) {
  double x = p[0];
  double y = p[1];
  return 2 * x + 4 * y;
}

// Source S = 0
double s_exact(const Amanzi::AmanziGeometry::Point& p) {
  double x = p[0];
  double y = p[1];
  return 0.0;
}

/* **************************************************************** */
TEST(MULTIPHASE_PRESSURE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase pressure equation, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  //std::string xmlFileName = "test/test_pressure_pk.xml";
  std::string xmlFileName = "test/test_pressure_spe10.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  std::string file_name = plist.sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  /*
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 5, 5, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  */
  RCP<const Mesh> mesh = meshfactory(file_name, gm);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(*state_list));
  //RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  std::string passwd("state");

  std::vector<Teuchos::RCP<Visualization> > visualization_;
  for (State::mesh_iterator mesh=S->mesh_begin();
       mesh!=S->mesh_end(); ++mesh) {
    if (parameter_list->isSublist("Visualization Data")) {
      Teuchos::ParameterList& vis_plist = parameter_list->sublist("Visualization Data");
      Teuchos::RCP<Visualization> vis =
        Teuchos::rcp(new Visualization(vis_plist, &comm));
      vis->set_mesh(mesh->second.first);
      vis->CreateFiles();
      visualization_.push_back(vis);
    }
  }
  
  const Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);
  Teuchos::RCP<CompositeVector> solution = Teuchos::rcp(new CompositeVector(*cvs));
  solution->PutScalar(0.0);

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  //soln->SetData(solution);

  // Create pressure PK
  Multiphase::Pressure_PK* PPK = new Multiphase::Pressure_PK(*pk_tree_list, parameter_list, S, soln);
  PPK->Setup();
  S->Setup();
  S->InitializeFields();

  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S.ptr());
  }

  // Init the PK and solve for pressure
  PPK->Initialize();
  PPK->CommitStep(0.0, 0.0);
  S->CheckAllFieldsInitialized();
  PPK->AdvanceStep(0.0, 1.0, false);
  PPK->CommitStep(0.0, 1.0);
  S->advance_time(1.0);
  S->advance_cycle();

  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S.ptr());
  }

  /*
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
  for (int c = 0; c < p.MyLength(); c++) 
  {
    const Point& xc = mesh->cell_centroid(c);
    CHECK_CLOSE(u_exact(xc), p[0][c], 1e-5);
  }

  // get the data for viz
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_pressure_pk.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }
  */
  delete PPK;
}
