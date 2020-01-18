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

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Reduced2p2c_PK.hh"
#include "OperatorDefs.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"


/* **************************************************************** */
TEST(MULTIPHASE_REDUCED_2P2C) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase component pk, uniform rectangular mesh" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/multiphase_reduced_2p2c.xml";
  auto plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  Teuchos::ParameterList pk_tree_list = plist->sublist("PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(plist, "state");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(plist, "cycle driver");
  Teuchos::RCP<Teuchos::ParameterList> MPMC_specs = Teuchos::sublist(plist, "MPMC Specs");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  // RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 200, 10);
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 50, 5);

  // create a simple state populate it
  RCP<State> S = rcp(new State(*state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  std::vector<Teuchos::RCP<Visualization> > visualization_;
  for (State::mesh_iterator mesh=S->mesh_begin(); mesh!=S->mesh_end(); ++mesh) {
    if (plist->isSublist("visualization data")) {
      Teuchos::ParameterList& vis_plist = plist->sublist("visualization data");
      Teuchos::RCP<Visualization> vis = Teuchos::rcp(new Visualization(vis_plist));
      vis->set_mesh(mesh->second.first);
      vis->CreateFiles();
      visualization_.push_back(vis);
    }
  }

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("Multiphase_PK", *plist));


  double dT = time_list->get<double>("time step");
  double t_final = time_list->get<double>("end time");

  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true)->SetOwned(false);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);

  Teuchos::RCP<CompositeVector> saturation = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> rhl = Teuchos::rcp(new CompositeVector(*cvs));
  saturation->PutScalar(1.0);
  pressure->PutScalar(1.0e6);
  rhl->PutScalar(0.0);

  Teuchos::RCP<CompositeVector> saturation_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> rhl_old = Teuchos::rcp(new CompositeVector(*cvs));

  Teuchos::RCP<TreeVector> sat_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> rhl_tree = Teuchos::rcp(new TreeVector());
  sat_tree->SetData(saturation);
  pres_tree->SetData(pressure);
  rhl_tree->SetData(rhl);

  Teuchos::RCP<TreeVector> sat_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> rhl_tree_old = Teuchos::rcp(new TreeVector());
  sat_tree_old->SetData(saturation_old);
  pres_tree_old->SetData(pressure_old);
  rhl_tree_old->SetData(rhl_old);

  Teuchos::RCP<CompositeVector> res1 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res2 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res3 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<TreeVector> f = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f3 = Teuchos::rcp(new TreeVector());

  f1->SetData(res1);
  f2->SetData(res2);
  f3->SetData(res3);
  f->PushBack(f1);
  f->PushBack(f2);
  f->PushBack(f3); 

  // create the global solution vector
  auto soln = Teuchos::rcp(new TreeVector());
  auto u_new = Teuchos::rcp(new TreeVector());
  u_new->PushBack(pres_tree);
  u_new->PushBack(sat_tree);
  u_new->PushBack(rhl_tree);

  auto u_old = Teuchos::rcp(new TreeVector());
  u_old->PushBack(pres_tree_old);
  u_old->PushBack(sat_tree_old);
  u_old->PushBack(rhl_tree_old);

  // create PK and initialize it
  auto MPMC = Teuchos::rcp(new Multiphase::Reduced2p2c_PK(pk_tree_list, plist, S, soln));
  MPMC->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  MPMC->Initialize(S.ptr());

  bool failed = true;
  double t_sim = 0.0;
  // double dT = 5000.0*365.0*24.0*3600.0;
  // double t_final = 10*dT;
  S->WriteStatistics(vo);

  while (!(t_sim > t_final || abs(t_sim - t_final) < 1e-8)) {
    if (t_sim + dT - t_final > 1e-12) dT = t_final - t_sim;
    failed = MPMC->AdvanceStep(t_sim, t_sim + dT, false);
    if (failed) {
      dT = MPMC->get_dt();
    } else {
      t_sim += dT;
      MPMC->CommitStep(t_sim, t_sim + dT, S); 
      S->advance_time(dT);
      S->advance_cycle();
      if (MyPID == 0) {
        std::cout << "State time: " << S->time() << "; State cycle: " << S->cycle() << endl;
      }
      for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin(); vis!=visualization_.end(); ++vis) {
        WriteVis((*vis).ptr(), S.ptr());
      } 
      dT = MPMC->get_dt();
    } 
  }
  S->WriteStatistics(vo);

  // MPMC->Functional(0.0, 0.1, u_new, u_new, f);
  // MPMC->UpdatePreconditioner(0.0, u_new, 0.1);

  // visualization
  std::string passwd("state");
  Epetra_MultiVector& pw = *S->GetFieldData("pressure_w", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& pn = *S->GetFieldData("pressure_n", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& sw = *S->GetFieldData("saturation_w", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& sn = *S->GetFieldData("saturation_n", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& hydrogen_density = *soln->SubVector(2)->Data()->ViewComponent("cell", false);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"multiphase_reduced_2p2c.gmv");
    GMV::start_data();
    GMV::write_cell_data(pw, 0, "pressure_w");
    GMV::write_cell_data(pn, 0, "pressure_n");
    GMV::write_cell_data(sw, 0, "saturation_w");
    GMV::write_cell_data(sn, 0, "saturation_n");
    GMV::write_cell_data(hydrogen_density, 0, "hydrogen_density");
    GMV::close_data_file();
  }
}
