/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Default.hh"
#include "PK_MixinLeaf.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinImplicitSubcycled.hh"
#include "PK_Adaptors.hh"

#include "pks_test.hh"

using namespace Amanzi;


//
// Creates an explicit PK
// ============================================================================
std::pair<Teuchos::RCP<PK>,Teuchos::RCP<State> >
createExplicit(const std::string& eqn_name, const std::string& ti_name) {
  std::string pk_name = eqn_name + ", " + ti_name;
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode_forward_euler.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));
  
  // intentionally leaks memory
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, comm));

  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  // make a 1x1 'mesh' for ODEs
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1, gm);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<PK> pk;
  if (eqn_name == "A") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorA> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "B") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "C") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> >(pk_tree, global_list, S, soln));
  } else {
    ASSERT(false);
  }

  return std::make_pair(pk, S);
}


//
// Creates an implicit PK
// ============================================================================
std::pair<Teuchos::RCP<PK>,Teuchos::RCP<State> >
createImplicit(const std::string& eqn_name) {
  std::string pk_name = eqn_name + ", backward euler";
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode_forward_euler.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));
  
  // intentionally leaks memory
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, comm));

  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  // make a 1x1 'mesh' for ODE
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1, gm);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<PK> pk;
  if (eqn_name == "A") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorA> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "B") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "C") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> >(pk_tree, global_list, S, soln));
  } else {
    ASSERT(false);
  }

  return std::make_pair(pk, S);
}


//
// Creates an implicit PK that always takes the recommended step size
// ============================================================================
std::pair<Teuchos::RCP<PK>,Teuchos::RCP<State> >
createImplicitSubcyled(const std::string& eqn_name) {
  std::string pk_name = eqn_name + ", backward euler subcycled";
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode_forward_euler.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));
  
  // intentionally leaks memory
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, comm));

  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  // make a 1x1 'mesh' for ODE
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1, gm);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<PK> pk;
  if (eqn_name == "A") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorA> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "B") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "C") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> >(pk_tree, global_list, S, soln));
  } else {
    ASSERT(false);
  }

  return std::make_pair(pk, S);
}


//
// Helper function that mocks a coordinator to run the test.
// ================================================================================
std::pair<double,double> run_test(const Teuchos::RCP<PK>& pk,
                                  const Teuchos::RCP<State>& S) {
  pk->Setup();
  S->Setup();
  pk->Initialize();
  S->Initialize();

  int nsteps_good = 0;
  int nsteps_bad = 0;

  double t_final = 1.0;
  bool done = false;
  while (!done) {
    double dt = pk->get_dt();
    S->set_time("next", S->time()+dt);
    S->set_cycle("next", S->cycle()+1);
    bool fail = pk->AdvanceStep("", "next");

    if (fail) {
      pk->FailStep("", "next");
      nsteps_bad++;
    } else {
      pk->CommitStep("", "next");
      S->set_time("", S->time("next"));
      S->set_cycle("", S->cycle("next"));
      nsteps_good++;
    }

    done = S->time("") >= t_final-0.0001;
  }
  return std::make_pair(nsteps_good, nsteps_bad);
}



SUITE(PK_ODE_FORWARD_EULER) {

  // Forward Euler tests with each of 3 PKs
  TEST(A_FORWARD_EULER) {
    auto pkS = createExplicit("A", "forward euler");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(2.0, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_FORWARD_EULER) {
    auto pkS = createExplicit("B", "forward euler");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(std::exp(1), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.59374, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_FORWARD_EULER) {
    auto pkS = createExplicit("C", "forward euler");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(std::exp(1), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.4);
    CHECK_CLOSE(2.33463, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }


  // Runge Kutta (multistage explicit) tests with each of 3 PKs
  TEST(A_RK4) {
    auto pkS = createExplicit("A", "RK4");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(2, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_RK4) {
    auto pkS = createExplicit("B", "RK4");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(std::exp(1), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_RK4) {
    auto pkS = createExplicit("C", "RK4");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(std::exp(1), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  
  // Implicit (single stage)
  TEST(A_BACKWARD_EULER) {
    auto pkS = createImplicit("A");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(2, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_BACKWARD_EULER) {
    auto pkS = createImplicit("B");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(std::exp(1.0), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.867971990790009, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER) {
    auto pkS = createImplicit("C");
    auto nsteps = run_test(pkS.first, pkS.second);
    CHECK_CLOSE(std::exp(1.0), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.6);
    CHECK_CLOSE(3.27476584420779, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }


  // Implicit (potentially multistage)
  TEST(A_BACKWARD_EULER_SUBCYCLED) {
    auto pkS = createImplicitSubcyled("A");
    auto nsteps = run_test(pkS.first, pkS.second);

    CHECK_CLOSE(2, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_BACKWARD_EULER_SUBCYCLED) {
    auto pkS = createImplicitSubcyled("B");
    auto nsteps = run_test(pkS.first, pkS.second);

    CHECK_CLOSE(std::exp(1.0), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.867971990790009, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER_SUBCYCLED) {
    auto pkS = createImplicitSubcyled("C");
    auto nsteps = run_test(pkS.first, pkS.second);

    CHECK_CLOSE(std::exp(1.0), (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.6);
    CHECK_CLOSE(3.27476584420779, (*pkS.second->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
  
}

