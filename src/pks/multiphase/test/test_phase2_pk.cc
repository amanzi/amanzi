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
#include "Phase2_PK.hh"
#include "OperatorDefs.hh"

double p_exact(const Amanzi::AmanziGeometry::Point& xc) {
  return -2.0*xc[0] + 2.0;
}

/* **************************************************************** */
TEST(MULTIPHASE_PHASE2) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase phase2, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  //std::string xmlFileName = "test/test_pressure_pk.xml";
  std::string xmlFileName = "test/test_phase2_pk.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");

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
  //cvs->AddComponent("face", AmanziMesh::FACE, 1);
  Teuchos::RCP<CompositeVector> saturation = Teuchos::rcp(new CompositeVector(*cvs));
  saturation->PutScalar(0.2);
  Teuchos::RCP<CompositeVector> pressure = Teuchos::rcp(new CompositeVector(*cvs));
  pressure->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> saturation_old = Teuchos::rcp(new CompositeVector(*cvs));
  saturation_old->PutScalar(0.2);
  Teuchos::RCP<CompositeVector> pressure_old = Teuchos::rcp(new CompositeVector(*cvs));
  pressure_old->PutScalar(0.0);
  Teuchos::RCP<CompositeVector> res = Teuchos::rcp(new CompositeVector(*cvs));
  res->PutScalar(0.0);
  /* 
  Epetra_MultiVector& saturation_cell = *saturation->ViewComponent("cell");
  Epetra_MultiVector& pressure_cell = *pressure->ViewComponent("cell");
  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point xc = mesh->cell_centroid(c);
    if (xc[0] < 0.35) {
      saturation_cell[0][c] = 1;
    }
    else {
      saturation_cell[0][c] = 1;
    }
    //pressure_cell[0][c] = p_exact(xc);
  }
  */

  Teuchos::RCP<TreeVector> sat_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree = Teuchos::rcp(new TreeVector());
  sat_tree->SetData(saturation);
  pres_tree->SetData(pressure);
  Teuchos::RCP<TreeVector> sat_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree_old = Teuchos::rcp(new TreeVector());
  sat_tree_old->SetData(saturation_old);
  pres_tree_old->SetData(pressure_old);
  Teuchos::RCP<TreeVector> f = Teuchos::rcp(new TreeVector());
  f->SetData(res);

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> u_new = Teuchos::rcp(new TreeVector());
  u_new->PushBack(pres_tree);
  u_new->PushBack(sat_tree);
  Teuchos::RCP<TreeVector> u_old = Teuchos::rcp(new TreeVector());
  u_old->PushBack(pres_tree_old);
  u_old->PushBack(sat_tree_old);
  //soln->SetData(solution);

  // Create pressure PK
  Multiphase::Phase2_PK* PPK = new Multiphase::Phase2_PK(*pk_tree_list, parameter_list, S, u_new);
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  PPK->Initialize();
  /*
  Epetra_MultiVector& phase_velocity = *S->GetFieldData("phase2_velocity", passwd)->ViewComponent("face", false);
  Point vel(1.0, 0.0);
  for (int f = 0; f < phase_velocity.MyLength(); f++) {
    phase_velocity[0][f] = vel * mesh->face_normal(f);
  }
  */
 
  PPK->UpdatePreconditioner(0.0, u_new, 0.1);
  PPK->Functional(0.0, 0.1, u_old, u_new, f);
  /*
  PPK->CommitStep(0.0, 0.0);
  S->CheckAllFieldsInitialized();
  PPK->AdvanceStep(0.0, 1.0);
  //PPK->CommitState(S.ptr());
  PPK->CommitStep(0.0, 1.0);

  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
  /*
  for (int c = 0; c < p.MyLength(); c++) 
  {
    const Point& xc = mesh->cell_centroid(c);
    CHECK_CLOSE(u_exact(xc), p[0][c], 1e-5);
  }
  */

  // get the data for viz
  Epetra_MultiVector& p = *S->GetFieldData("phase2_pressure", passwd)->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_phase2_pk.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }
  delete PPK;
}
