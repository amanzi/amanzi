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
#include "PK_MixinExplicitSubcycled.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinImplicitSubcycled.hh"
#include "PK_Adaptors.hh"

#include "pks_test.hh"

using namespace Amanzi;


struct Run {
  Run(const Teuchos::RCP<State>& S_,
      const Teuchos::RCP<PK>& pk_,
      const Teuchos::RCP<TreeVector>& soln_)
      : S(S_),
        pk(pk_),
        soln(soln_) {}

  Teuchos::RCP<State> S;
  Teuchos::RCP<PK> pk;
  Teuchos::RCP<TreeVector> soln;
};

//
// Creates an explicit PK
// ============================================================================
std::unique_ptr<Run>
createExplicit(const std::string& eqn_name, const std::string& ti_name) {
  std::string pk_name = eqn_name + ", " + ti_name;
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
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

  return std::make_unique<Run>(S,pk,soln);
}


//
// Creates an explicit subcycled PK
// ============================================================================
std::unique_ptr<Run>
createExplicitSubcycled(const std::string& eqn_name, const std::string& ti_name) {
  std::string pk_name = eqn_name + ", " + ti_name + " subcycled";
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
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
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorA> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "B") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> >(pk_tree, global_list, S, soln));
  } else if (eqn_name == "C") {
    pk = Teuchos::rcp(new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> >(pk_tree, global_list, S, soln));
  } else {
    ASSERT(false);
  }

  return std::make_unique<Run>(S,pk,soln);
}


//
// Creates an implicit PK
// ============================================================================
std::unique_ptr<Run>
createImplicit(const std::string& eqn_name) {
  std::string pk_name = eqn_name + ", backward euler";
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
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

  return std::make_unique<Run>(S,pk,soln);
}



//
// Creates a subcycled implicit PK.  These don't fail, they simply subcycle as
// needed to match the full timestep.
// ============================================================================
std::unique_ptr<Run>
createImplicitSubcycled(const std::string& eqn_name) {
  std::string pk_name = eqn_name + ", backward euler subcycled";
  std::cout << "Test: " << pk_name << std::endl;
  
  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
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

  return std::make_unique<Run>(S,pk,soln);
}



//
// Helper function that mocks a coordinator to run the test.
// ================================================================================
std::pair<double,double> run_test(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<PK>& pk,
        const Teuchos::RCP<TreeVector>& soln_next) {

  auto soln_old = Teuchos::rcp(new TreeVector(*soln_next, INIT_MODE_NOALLOC));
  pk->SolutionToState(*soln_old, "", "");
  pk->SolutionToState(*soln_next, "next", "");

  pk->Setup(*soln_next);
  S->Setup();

  pk->StateToSolution(*soln_old, "", "");
  pk->StateToSolution(*soln_next, "next", "");

  pk->Initialize();
  S->Initialize();

  *soln_next = *soln_old;

  int nsteps_good = 0;
  int nsteps_bad = 0;

  double t_final = 1.0;
  bool done = false;
  while (!done) {
    double dt = std::min(pk->get_dt(), 1.0);
    S->set_time("next", S->time()+dt);
    S->set_cycle("next", S->cycle()+1);

    bool fail = pk->AdvanceStep("", soln_old, "next", soln_next);

    if (fail) {
      pk->FailStep("", soln_old, "next", soln_next);
      nsteps_bad++;
    } else {
      pk->CommitStep("", soln_old, "next", soln_next);
      S->set_time("", S->time("next"));
      S->set_cycle("", S->cycle("next"));
      nsteps_good++;
    }

    done = S->time("") >= t_final-0.0001;
  }
  return std::make_pair(nsteps_good, nsteps_bad);
}



SUITE(PKS_ODE) {

  // Forward Euler tests with each of 3 PKs
  TEST(A_FORWARD_EULER) {
    auto run = createExplicit("A", "forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(2.0, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-10);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_FORWARD_EULER) {
    auto run = createExplicit("B", "forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.59374, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_FORWARD_EULER) {
    auto run = createExplicit("C", "forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.4);
    CHECK_CLOSE(2.33463, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }


  // Runge Kutta (multistage explicit) tests with each of 3 PKs
  TEST(A_RK4) {
    auto run = createExplicit("A", "RK4");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(2, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_RK4) {
    auto run = createExplicit("B", "RK4");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_RK4) {
    auto run = createExplicit("C", "RK4");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  // Subcycled Forward Euler tests with each of 3 PKs.  Note that this answer
  // is identical to the Forward Euler -- it takes 10 subcycled steps for an
  // outer step that is 10 times as big.
  TEST(A_FORWARD_EULER_SUBCYCLED) {
    auto run = createExplicitSubcycled("A", "forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(2.0, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-10);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_FORWARD_EULER_SUBCYCLED) {
    auto run = createExplicitSubcycled("B", "forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.59374, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_FORWARD_EULER_SUBCYCLED) {
    auto run = createExplicitSubcycled("C", "forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.4);
    CHECK_CLOSE(2.33463, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
  
  TEST(A_RK4_SUBCYCLED) {
    auto run = createExplicitSubcycled("A", "RK4");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(2.0, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-10);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }


  TEST(B_RK4_SUBCYCLED) {
    auto run = createExplicitSubcycled("B", "RK4");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_RK4_SUBCYCLED) {
    auto run = createExplicitSubcycled("C", "RK4");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
  

  // Implicit (single stage)
  TEST(A_BACKWARD_EULER) {
    auto run = createImplicit("A");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(2, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_BACKWARD_EULER) {
    auto run = createImplicit("B");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.867971990790009, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER) {
    auto run = createImplicit("C");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.6);
    CHECK_CLOSE(3.27476584420779, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  
  // Implicit (subcycled)
  TEST(A_BACKWARD_EULER_SUBCYCLED) {
    auto run = createImplicitSubcycled("A");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(2, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_BACKWARD_EULER_SUBCYCLED) {
    auto run = createImplicitSubcycled("B");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.86608, (*run->S->Get<CompositeVector>("primary")
            .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER_SUBCYCLED) {
    auto run = createImplicitSubcycled("C");
    auto nsteps = run_test(run->S, run->pk, run->soln);
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 0.3);
    CHECK_CLOSE(3.01043, (*run->S->Get<CompositeVector>("primary")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
  
}

