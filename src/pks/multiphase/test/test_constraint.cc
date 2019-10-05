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
#include "Phase_Constraint_PK.hh"
#include "OperatorDefs.hh"


/* **************************************************************** */
TEST(MULTIPHASE_COMP) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase component pk, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  //std::string xmlFileName = "test/test_pressure_pk.xml";
  std::string xmlFileName = "test/test_comp_pk.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  Teuchos::RCP<Teuchos::ParameterList> comp_list = Teuchos::rcp(new Teuchos::ParameterList(plist.sublist("MPMC Specs")));

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 0.1, 5, 1, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  /* create a simple state populate it */
  //RCP<State> S = rcp(new State(*state_list));
  RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  std::string passwd("state");
  
  const Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);

  Teuchos::RCP<CompositeVector> saturation = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> fug_a = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> fug_w = Teuchos::rcp(new CompositeVector(*cvs));
  saturation->PutScalar(0.2);
  pressure->PutScalar(10.0);
  fug_a->PutScalar(6.0);
  fug_w->PutScalar(4.0);

  Teuchos::RCP<CompositeVector> saturation_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> fug_a_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> fug_w_old = Teuchos::rcp(new CompositeVector(*cvs));
  saturation_old->PutScalar(0.2);
  pressure_old->PutScalar(10.0);
  fug_a_old->PutScalar(8.0);
  fug_w_old->PutScalar(2.0);


  Teuchos::RCP<TreeVector> sat_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> fug_a_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> fug_w_tree = Teuchos::rcp(new TreeVector());
  sat_tree->SetData(saturation);
  pres_tree->SetData(pressure);
  fug_a_tree->SetData(fug_a);
  fug_w_tree->SetData(fug_w);

  Teuchos::RCP<TreeVector> sat_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> fug_a_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> fug_w_tree_old = Teuchos::rcp(new TreeVector());
  sat_tree_old->SetData(saturation_old);
  pres_tree_old->SetData(pressure_old);
  fug_a_tree_old->SetData(fug_a_old);
  fug_w_tree_old->SetData(fug_w_old);

  Teuchos::RCP<CompositeVector> res1 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res2 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res3 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res4 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f3 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f4 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f = Teuchos::rcp(new TreeVector());
  f1->SetData(res1);
  f2->SetData(res2);
  f3->SetData(res3);
  f4->SetData(res4);
  f->PushBack(f1);
  f->PushBack(f2);
  f->PushBack(f3); 
  f->PushBack(f4);

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> u_new = Teuchos::rcp(new TreeVector());
  u_new->PushBack(pres_tree);
  u_new->PushBack(sat_tree);
  u_new->PushBack(fug_w_tree);
  u_new->PushBack(fug_a_tree);
  Teuchos::RCP<TreeVector> u_old = Teuchos::rcp(new TreeVector());
  u_old->PushBack(pres_tree_old);
  u_old->PushBack(sat_tree_old);
  u_old->PushBack(fug_w_tree_old);
  u_old->PushBack(fug_a_tree_old);

  // Create pressure PK
  Multiphase::Phase_Constraint_PK* PC_PK = new Multiphase::Phase_Constraint_PK(*comp_list, S);
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  PC_PK->Initialize();
  PC_PK->Functional(0.0, 0.1, u_old, u_new, f); // phase constraints are time independent
  std::cout << "Functional f3: " << *f->SubVector(2)->Data()->ViewComponent("cell") << "\n";
  std::cout << "Functional f4: " << *f->SubVector(3)->Data()->ViewComponent("cell") << "\n"; 
  PC_PK->UpdatePreconditioner(0.0, u_old, 0.1);

  // get the data for viz
  /*
  Epetra_MultiVector& p = *S->GetFieldData("phase_pressure", passwd)->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_phase1_pk.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }
  */
  delete PC_PK;
}
