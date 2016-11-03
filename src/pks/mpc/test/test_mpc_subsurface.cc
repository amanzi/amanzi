/*
  Testing of subsurface coupled thermal hydrology

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "MeshFactory.hh"
//#include "mpc_subsurface.hh"
#include "weak_mpc.hh"
#include "three_phase.hh"
#include "permafrost.hh"
#include "TreeOperator.hh"
#include "State.hh"

/* **************************************************************** */
TEST(MPC_SUBSURFACE_HYDROLOGY) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: hydrology" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/hydrology.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> global_list = Teuchos::rcpFromRef(plist);

  // create a mesh framework
  Teuchos::RCP<Teuchos::ParameterList> region_list = Teuchos::sublist(global_list, "Regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, *region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory(0.,0.,0., 1.,1.,10., 1,1,100, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(global_list, "state");

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_list, "PK");
  Teuchos::RCP<Teuchos::ParameterList> fe_list = Teuchos::sublist(state_list, "field evaluators");


  Teuchos::RCP<State> S = Teuchos::rcp(new State(*state_list));
  S->RegisterDomainMesh(mesh);

  Flow::Permafrost mpc(*fe_list, pk_list,  S, soln);
 
  mpc.Setup(S.ptr());

  Teuchos::ParameterList temp_eval_plist;
  temp_eval_plist.set("evaluator name", "temperature");
  Teuchos::RCP<PrimaryVariableFieldEvaluator> temp_pv = Teuchos::rcp(new PrimaryVariableFieldEvaluator(temp_eval_plist));
  S->SetFieldEvaluator("temperature", temp_pv);
  S->RequireField("temperature", "energy")->SetMesh(mesh)->SetGhosted()
      ->AddComponent("cell",AmanziMesh::CELL,1)->AddComponent("face",AmanziMesh::FACE,1);

  S->Setup();

  // initialize manually
  Teuchos::ParameterList pressure_list;
  pressure_list.set("restart file", "test/checkpoint-mpc_subsurface.h5");
  S->GetField("pressure", "flow")->Initialize(pressure_list);

  Teuchos::ParameterList temp_list;
  temp_list.set("restart file", "test/checkpoint-mpc_subsurface.h5");
  S->GetField("temperature", "energy")->Initialize(temp_list);
  
  mpc.Initialize(S.ptr());
  S->Initialize();
  S->CheckAllFieldsInitialized();

  Teuchos::RCP<State> S_old = Teuchos::rcp(new State(*S));
  *S_old = *S;
  mpc.set_states(S_old, S_old, S);

  // check the preconditioner
  // -- ApplyInverse(Apply()) is the identity
  Teuchos::RCP<TreeVector> du = Teuchos::rcp(new TreeVector(*soln));
  Teuchos::RCP<TreeVector> du_test = Teuchos::rcp(new TreeVector(*soln));
  Teuchos::RCP<TreeVector> r = Teuchos::rcp(new TreeVector(*soln));
  //  du->Data()->Random();
  du->PutScalar(1.);

  mpc.UpdatePreconditioner(0., soln, 864.);
  mpc.preconditioner()->Apply(*du->Data(), *r->Data());
  // std::cout << "r_flow" << std::endl;
  // r->Data()->Print(std::cout);
  mpc.ApplyPreconditioner(r, du_test);
  du_test->Update(-1, *du, 1);
  double enorm;
  du_test->NormInf(&enorm);
  CHECK_CLOSE(0., enorm, 1.e-10);

}


/* **************************************************************** */
TEST(MPC_SUBSURFACE_ENERGY) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: energy" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/energy.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> global_list = Teuchos::rcpFromRef(plist);

  // create a mesh framework
  Teuchos::RCP<Teuchos::ParameterList> region_list = Teuchos::sublist(global_list, "Regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, *region_list, &comm));

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory(0.,0.,0., 1.,1.,10., 1,1,100, gm); 

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(global_list, "state");

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_list, "PK");
  Teuchos::RCP<Teuchos::ParameterList> fe_list = Teuchos::sublist(state_list, "field evaluators");


  Teuchos::RCP<State> S = Teuchos::rcp(new State(*state_list));
  S->RegisterDomainMesh(mesh);
  S->RequireScalar("atmospheric_pressure");


  Energy::ThreePhase mpc(*fe_list, pk_list, S, soln); 
  mpc.Setup(S.ptr());

  Teuchos::ParameterList pres_eval_plist;
  pres_eval_plist.set("evaluator name", "pressure");
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pres_pv = Teuchos::rcp(new PrimaryVariableFieldEvaluator(pres_eval_plist));
  S->SetFieldEvaluator("pressure", pres_pv);
  S->RequireField("pressure", "flow")->SetMesh(mesh)->SetGhosted()
      ->AddComponent("cell",AmanziMesh::CELL,1)->AddComponent("face",AmanziMesh::FACE,1);

  S->RequireField("mass_flux", "mass_flux")->SetMesh(mesh)->SetGhosted()
      ->AddComponent("face",AmanziMesh::FACE,1);
  
  S->Setup();

  // initialize manually
  Teuchos::ParameterList pressure_list;
  pressure_list.set("restart file", "test/checkpoint-mpc_subsurface.h5");
  S->GetField("pressure", "flow")->Initialize(pressure_list);

  Teuchos::ParameterList temp_list;
  temp_list.set("restart file", "test/checkpoint-mpc_subsurface.h5");
  S->GetField("temperature", "energy")->Initialize(temp_list);

  S->GetFieldData("mass_flux", "mass_flux")->PutScalar(0.);
  S->GetField("mass_flux", "mass_flux")->set_initialized();
  
  mpc.Initialize(S.ptr());
  S->Initialize();
  S->CheckAllFieldsInitialized();

  Teuchos::RCP<State> S_old = Teuchos::rcp(new State(*S));
  *S_old = *S;
  mpc.set_states(S_old, S_old, S);

  // check the preconditioner
  // -- ApplyInverse(Apply()) is the identity
  Teuchos::RCP<TreeVector> du = Teuchos::rcp(new TreeVector(*soln));
  Teuchos::RCP<TreeVector> du_test = Teuchos::rcp(new TreeVector(*soln));
  Teuchos::RCP<TreeVector> r = Teuchos::rcp(new TreeVector(*soln));
  //  du->Data()->Random();
  du->PutScalar(1.);

  mpc.UpdatePreconditioner(0., soln, 864.);
  mpc.preconditioner()->Apply(*du->Data(), *r->Data());

  // std::cout << "r_energy" << std::endl;
  // r->Data()->Print(std::cout);

  mpc.ApplyPreconditioner(r, du_test);
  du_test->Update(-1, *du, 1);
  double enorm;
  du_test->NormInf(&enorm);
  CHECK_CLOSE(0., enorm, 1.e-10);

}



// // // /* **************************************************************** */
// // TEST(MPC_SUBSURFACE_THERMAL_HYDROLOGY) {
// //   using namespace Amanzi;
// //   using namespace Amanzi::AmanziGeometry;
// //   using namespace Amanzi::AmanziMesh;

// //   Epetra_MpiComm comm(MPI_COMM_WORLD);
// //   int MyPID = comm.MyPID();

// //   if (MyPID == 0) std::cout << "Test: coupled thermal hydrology" << std::endl;

// //   /* read parameter list */
// //   std::string xmlFileName = "test/mpc_subsurface.xml";
// //   Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
// //   Teuchos::ParameterList plist = xmlreader.getParameters();
// //   Teuchos::RCP<Teuchos::ParameterList> global_list = Teuchos::rcpFromRef(plist);

// //   // create a mesh framework
// //   Teuchos::RCP<Teuchos::ParameterList> region_list = Teuchos::sublist(global_list, "Regions");
// //   Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, *region_list, &comm));

// //   FrameworkPreference pref;
// //   pref.clear();
// //   pref.push_back(MSTK);
// //   pref.push_back(STKMESH);

// //   MeshFactory meshfactory(&comm);
// //   meshfactory.preference(pref);
// //   Teuchos::RCP<Mesh> mesh = meshfactory(0.,0.,0., 1.,1.,10., 1,1,100, gm); 

// //   /* create a simple state and populate it */
// //   Amanzi::VerboseObject::hide_line_prefix = true;

// //   Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(global_list, "state");

// //   Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
// //   Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_list, "PK");
// //   Teuchos::RCP<Teuchos::ParameterList> fe_list = Teuchos::sublist(state_list, "field evaluators");
// //   MPCSubsurface mpc(pk_list, *fe_list, soln);

// //   Teuchos::RCP<State> S = Teuchos::rcp(new State(*state_list));
// //   S->RegisterDomainMesh(mesh);
 
// //   mpc.setup(S.ptr());
// //   S->Setup();

// //   // initialize manually
// //   Teuchos::ParameterList pressure_list;
// //   pressure_list.set("restart file", "test/checkpoint-mpc_subsurface.h5");
// //   S->GetField("pressure", "flow")->Initialize(pressure_list);

// //   Teuchos::ParameterList temp_list;
// //   temp_list.set("restart file", "test/checkpoint-mpc_subsurface.h5");
// //   S->GetField("temperature", "energy")->Initialize(temp_list);  
  
// //   mpc.initialize(S.ptr());
// //   S->Initialize();
// //   S->CheckAllFieldsInitialized();

// //   Teuchos::RCP<State> S_old = Teuchos::rcp(new State(*S));
// //   *S_old = *S;
// //   mpc.set_states(S_old, S_old, S);

// //   // check the preconditioner
// //   // -- ApplyInverse(Apply()) is the identity
// //   Teuchos::RCP<TreeVector> du = Teuchos::rcp(new TreeVector(*soln));
// //   Teuchos::RCP<TreeVector> du_test = Teuchos::rcp(new TreeVector(*soln));
// //   Teuchos::RCP<TreeVector> r = Teuchos::rcp(new TreeVector(*soln));
// //   // du->SubVector(0)->Data()->Random();
// //   // du->SubVector(1)->Data()->Random();
// //   du_test->PutScalar(0.);
// //   //  du->PutScalar(1.);
// //   du->SubVector(0)->PutScalar(1.);
// //   du->SubVector(1)->PutScalar(1.);

// //   mpc.UpdatePreconditioner(0., soln, 864.);
// //   mpc.preconditioner()->Apply(*du, *r);
// //   // std::cout << "r_coupled:" << std::endl;
// //   // r->Print(std::cout);
// //   mpc.ApplyPreconditioner(r, du_test);
// //   // std::cout << "inverted solution (should be 1):" << std::endl;
// //   // du_test->Print(std::cout);

// //   du_test->Update(-1, *du, 1);
// //   // std::cout << "inverted solution - soln (should be 0):" << std::endl;
// //   // du_test->Print(std::cout);
// //   double enorm_e, enorm_f;
// //   du_test->SubVector(0)->NormInf(&enorm_f);
// //   du_test->SubVector(1)->NormInf(&enorm_e);
// //   double enorm = std::max(enorm_e, enorm_f);

// //   CHECK_CLOSE(0., enorm, 1.e-8);

// // }
