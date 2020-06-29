/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Adaptors.hh"
#include "PK_Default.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"

using namespace Amanzi;

template <class Base_t>
class PK_ODE_Explicit : public Base_t {
 public:
  PK_ODE_Explicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution)
    : Base_t(pk_tree, global_plist, S, solution)
  {}

  void Setup()
  {
    Base_t::Setup();

    this->S_
      ->template Require<CompositeVector, CompositeVectorSpace>(this->key_)
      .SetComponent("cell", AmanziMesh::CELL, 1)
      ->SetGhosted(false);
  }

  void Initialize()
  {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>("primary", "", "primary")
      .PutScalar(1.);
    this->S_->GetRecordW("primary", "primary").set_initialized();
  }

  void FunctionalTimeDerivative(double t, const TreeVector& u, TreeVector& f)
  {
    f = u;
    std::cout << "  At time t = " << t
              << ": u = " << (*u.Data()->ViewComponent("cell", false))[0][0]
              << ", du/dt = " << (*f.Data()->ViewComponent("cell", false))[0][0]
              << std::endl;
  }
};

Teuchos::RCP<PK>
create(const Teuchos::RCP<State>& S,
       const Teuchos::RCP<Teuchos::ParameterList>& global_list)
{
  Teuchos::RCP<Teuchos::ParameterList> pk_tree =
    Teuchos::rcp(new Teuchos::ParameterList("my pk"));

  // intentionally leaks memory
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // create mesh
  Teuchos::ParameterList regions_list;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, comm));

  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
    meshfactory(0.0, -2.0, 1.0, 0.0, 18, 18, gm);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  return Teuchos::rcp(
    new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default>>>>(
      pk_tree, global_list, S, soln));
}

SUITE(PK_ODE_EXPLICIT2)
{
  TEST(Construct)
  {
    std::cout << "PK_ODE_Explicit2: Construction" << std::endl;

    Teuchos::RCP<Teuchos::ParameterList> global_list =
      Teuchos::rcp(new Teuchos::ParameterList("main"));
    auto sublist =
      Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
    sublist->set<std::string>("domain name", "domain");
    sublist->set<std::string>("primary variable key", "primary");

    auto S = Teuchos::rcp(new State());
    auto pk = create(S, global_list);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();
  }

  TEST(Advance)
  {
    std::cout << "PK_ODE_Explicit2: Advance 1 step" << std::endl;

    Teuchos::RCP<Teuchos::ParameterList> global_list =
      Teuchos::rcp(new Teuchos::ParameterList("main"));
    auto sublist =
      Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
    sublist->set<std::string>("domain name", "domain");
    sublist->set<std::string>("primary variable key", "primary");

    auto S = Teuchos::rcp(new State());
    auto pk = create(S, global_list);

    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    bool fail = pk->AdvanceStep("", "next");

    CHECK(!fail);
    CHECK_CLOSE(1.0,
                (*S->Get<CompositeVector>("primary", "")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
    CHECK_CLOSE(2.0,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
  }

  TEST(Advance2)
  {
    std::cout << "PK_ODE_Explicit2: Advance to finish" << std::endl;
    Teuchos::RCP<Teuchos::ParameterList> global_list =
      Teuchos::rcp(new Teuchos::ParameterList("main"));
    auto sublist =
      Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
    sublist->set<std::string>("domain name", "domain");
    sublist->set<std::string>("primary variable key", "primary");

    auto S = Teuchos::rcp(new State());
    auto pk = create(S, global_list);

    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double t_final = 1.0;
    while (S->time() < t_final) {
      double dt = pk->get_dt();
      S->set_time("next", S->time() + dt);
      S->set_cycle("next", S->cycle() + 1);
      bool fail = pk->AdvanceStep("", "next");

      if (fail) {
        pk->FailStep("", "next");
      } else {
        pk->CommitStep("", "next");
        S->set_time("", S->time("next"));
        S->set_cycle("", S->cycle("next"));
      }
    }

    CHECK_CLOSE(2.0,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
  }

  TEST(Advance_RK2)
  {
    std::cout << "PK_ODE_Explicit2: Advance to finish with RK midpoint"
              << std::endl;
    Teuchos::RCP<Teuchos::ParameterList> global_list =
      Teuchos::rcp(new Teuchos::ParameterList("main"));
    auto sublist =
      Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
    sublist->set<std::string>("domain name", "domain");
    sublist->set<std::string>("primary variable key", "primary");

    auto ti_sublist = Teuchos::sublist(sublist, "time integrator");
    ti_sublist->set<std::string>("RK method", "midpoint");

    auto S = Teuchos::rcp(new State());
    auto pk = create(S, global_list);

    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double t_final = 1.0;
    while (S->time() < t_final) {
      double dt = pk->get_dt();
      S->set_time("next", S->time() + dt);
      S->set_cycle("next", S->cycle() + 1);
      bool fail = pk->AdvanceStep("", "next");

      if (fail) {
        pk->FailStep("", "next");
      } else {
        pk->CommitStep("", "next");
        S->set_time("", S->time("next"));
        S->set_cycle("", S->cycle("next"));
      }
    }

    CHECK_CLOSE(2.5,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
  }

  TEST(Advance_RK4)
  {
    std::cout << "PK_ODE_Explicit2: Advance to finish with RK 4th order"
              << std::endl;
    Teuchos::RCP<Teuchos::ParameterList> global_list =
      Teuchos::rcp(new Teuchos::ParameterList("main"));
    auto sublist =
      Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
    sublist->set<std::string>("domain name", "domain");
    sublist->set<std::string>("primary variable key", "primary");
    sublist->set<double>("initial time step", 0.1);

    auto ti_sublist = Teuchos::sublist(sublist, "time integrator");
    ti_sublist->set<std::string>("RK method", "runge kutta 4th order");

    auto S = Teuchos::rcp(new State());
    auto pk = create(S, global_list);

    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double t_final = 1.0;
    bool done = false;
    while (!done) {
      double dt = pk->get_dt();
      S->set_time("next", S->time() + dt);
      S->set_cycle("next", S->cycle() + 1);
      bool fail = pk->AdvanceStep("", "next");

      if (fail) {
        pk->FailStep("", "next");
      } else {
        pk->CommitStep("", "next");
        S->set_time("", S->time("next"));
        S->set_cycle("", S->cycle("next"));
      }

      done = S->time("") >= t_final - 0.0001;
    }

    // solution, y = exp(t)|_{t=1} = e
    CHECK_CLOSE(2.7182818284590451,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-5);
  }
}
