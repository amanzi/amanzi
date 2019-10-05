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
#include "MPMC_PK.hh"
#include "OperatorDefs.hh"


/* **************************************************************** */
TEST(MULTIPHASE_COMP) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase component pk, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  //std::string xmlFileName = "test/test_pressure_pk.xml";
  std::string xmlFileName = "test/test_mpmc.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::ParameterList pk_tree_list = plist.sublist("PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector();
  std::vector<double> high_coor = mesh_list.get<Teuchos::Array<double> > ("Domain High Coordinate").toVector();
  std::vector<double> low_coor = mesh_list.get<Teuchos::Array<double> > ("Domain Low Coordinate").toVector();

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
  //RCP<State> S = rcp(new State());
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
  fug_a->PutScalar(8.0);
  fug_w->PutScalar(2.0);

  Teuchos::RCP<CompositeVector> saturation_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> fug_a_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> fug_w_old = Teuchos::rcp(new CompositeVector(*cvs));

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
  Teuchos::RCP<TreeVector> f = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f3 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f4 = Teuchos::rcp(new TreeVector());
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
  Multiphase::MPMC_PK* MPMC = new Multiphase::MPMC_PK(pk_tree_list, parameter_list, S, soln);
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  MPMC->Initialize();

  /*
  bool failed = true;
  double dT, T0, t_sim, t_fin;
  dT = 0.005;
  T0 = 0.0;
  t_sim = 0.0;

  failed = MPMC->AdvanceStep(T0, T0 + dT, false);
  if (!failed) MPMC->CommitStep(T0, T0 + dT);
  */

  bool failed = true;
  //double dT = 0.05;
  //double dT = 0.6/mesh_size[0];
  double dT = 0.01;
  double t_sim = 0.0;
  double t_final = 0.11;
  while (!(t_sim > t_final || abs(t_sim - t_final) < 1e-8)) {
    failed = MPMC->AdvanceStep(t_sim, t_sim + dT, false);
    if (failed) {
      dT /= 2.0;
    } else {
      t_sim += dT;
      MPMC->CommitStep(t_sim, t_sim + dT); 
    } 
  }
  std::cout << "t_sim: " << t_sim << "\n"; 

  //MPMC->Functional(0.0, 0.1, u_new, u_new, f);
  //MPMC->UpdatePreconditioner(0.0, u_new, 0.1);

  // get the data for viz
  Epetra_MultiVector& pw = *S->GetFieldData("pressure_w", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& pn = *S->GetFieldData("pressure_n", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& sw = *S->GetFieldData("saturation_w", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& sn = *S->GetFieldData("saturation_n", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& fuga1 = *soln->SubVector(2)->Data()->ViewComponent("cell", false);
  Epetra_MultiVector& fuga2 = *soln->SubVector(3)->Data()->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_mpmc.gmv");
    GMV::start_data();
    GMV::write_cell_data(pw, 0, "pressure_w");
    GMV::write_cell_data(pn, 0, "pressure_n");
    GMV::write_cell_data(sw, 0, "saturation_w");
    GMV::write_cell_data(sn, 0, "saturation_n");
    GMV::write_cell_data(fuga1, 0, "fugacity1");
    GMV::write_cell_data(fuga2, 0, "fugacity2");
    GMV::close_data_file();
  }

  delete MPMC;
}
