/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
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
#include "PK_MixinImplicit.hh"
#include "PK_MixinLeaf.hh"

using namespace Amanzi;

template <class Base_t>
class PK_ODE_Implicit : public Base_t {
 public:
  PK_ODE_Implicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
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
    this->S_->GetRecordW("primary", "", "primary").set_initialized();
  }

  void FunctionalResidual(double t, const TreeVector& u, TreeVector& f)
  {
    f = u;
  }

  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f)
  {
    double dt = t_new - t_old;
    FunctionalResidual(t_new, *u_new, *f);
    f->Update(1. / dt, *u_new, -1. / dt, *u_old, -1.);
    std::cout << "  At t = " << t_old << ": u_old = "
              << (*u_old->Data()->ViewComponent("cell", false))[0][0]
              << ", u_new = "
              << (*u_new->Data()->ViewComponent("cell", false))[0][0]
              << ", r = " << (*f->Data()->ViewComponent("cell", false))[0][0]
              << std::endl;
  }

  // applies preconditioner to u and returns the result in Pu
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                          Teuchos::RCP<TreeVector> Pu)
  {
    *Pu = *u;
    Pu->Scale(1.0 / h_);
    return 0;
  }

  // updates the preconditioner
  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
  {
    h_ = 1.0 / h - 1.0;
    if (std::abs(h_) < 1.e-6) h_ = 1.e-6;
  }

  // computes a norm on u-du and returns the result
  double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
  {
    double norm;
    du->Norm2(&norm);
    return norm;
  }

  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) { return true; }

  bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u)
  {
    return false;
  }

  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du)
  {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  void ChangedSolution() {}

 protected:
  double h_;
};

Teuchos::RCP<PK>
create(const Teuchos::RCP<State>& S)
{
  Teuchos::RCP<Teuchos::ParameterList> global_list =
    Teuchos::rcp(new Teuchos::ParameterList("main"));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree =
    Teuchos::rcp(new Teuchos::ParameterList("my pk"));
  auto sublist =
    Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), "my pk");
  sublist->set<std::string>("domain name", "domain");
  sublist->set<std::string>("primary variable key", "primary");
  sublist->set<double>("initial time step", 1.0);

  auto ti_list = Teuchos::sublist(sublist, "time integrator");
  ti_list->set<std::string>("solver type", "Newton");
  auto ti_s_list = Teuchos::sublist(ti_list, "Newton parameters");
  ti_s_list->set<double>("nonlinear tolerance", 1.e-12);
  ti_list->set<std::string>("timestep controller type", "standard");
  auto ti_c_list =
    Teuchos::sublist(ti_list, "timestep controller standard parameters");
  ti_c_list->set("min iterations", 1);
  ti_c_list->set("max iterations", 5);
  ti_c_list->set("time step reduction factor", 0.5);
  ti_c_list->set("time step increase factor", 2.0);
  ti_c_list->set("max time step", 2.0);
  ti_c_list->set("min time step", 0.25);

  // intentionally leaks memory
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_transient_2D.xml";

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
    new PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeaf<PK_Default>>>>(
      pk_tree, global_list, S, soln));
}

SUITE(PK_ODE_IMPLICIT_WITH_FAIL)
{
  TEST(Construct)
  {
    std::cout << "PK_ODE_Implicit2: Construct" << std::endl;
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();
  }

  TEST(Advance)
  {
    std::cout << "PK_ODE_Implicit2: Advance 1 step" << std::endl;
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    bool fail = pk->AdvanceStep("", "next");

    CHECK(fail);
  }

  TEST(Advance2)
  {
    std::cout << "PK_ODE_Implicit2: Advance to finish" << std::endl;
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
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

    CHECK_CLOSE(4.0,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
  }
}
