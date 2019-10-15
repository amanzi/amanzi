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
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"

using namespace Amanzi;

Teuchos::RCP<PK>
createForwardEuler(const Teuchos::RCP<State>& S, const std::string& eqn_name)
{
  std::string pk_name = eqn_name + ", forward euler" std::cout
                        << "Test: " << pk_name << std::endl;

  Teuchos::RCP<Teuchos::ParameterList> pk_tree =
    Teuchos::rcp(new Teuchos::ParameterList(pk_name));
  Teuchos::RCP<Teuchos::ParameterList> global_list =
    Teuchos::rcp(new Teuchos::ParameterList("main"));
  auto sublist =
    Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), pk_name);
  sublist->set<std::string>("domain name", "domain");
  sublist->set<std::string>("primary variable key", "primary");

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
  // make a 1x1 'mesh' for ODEs
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
    meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1, gm);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  if (eqn_name == "A") {
    return Teuchos::rcp(
      new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default>>,
                                     DudtEvaluatorA>>(
        pk_tree, global_list, S, soln));
  } else if (eqn_name == "B") {
    return Teuchos::rcp(
      new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default>>,
                                     DudtEvaluatorB>>(
        pk_tree, global_list, S, soln));
  } else if (eqn_name == "C") {
    return Teuchos::rcp(
      new PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default>>,
                                     DudtEvaluatorC>>(
        pk_tree, global_list, S, soln));
  } else {
    AMANZI_ASSERT(false);
  }
}

SUITE(PK_ODE_EXPLICIT)
{
  TEST(Construct)
  {
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();
  }

  TEST(Advance)
  {
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    double dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    pk->AdvanceStep("", "next");

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
    auto S = Teuchos::rcp(new State());
    auto pk = create(S);
    pk->Setup();
    S->Setup();
    pk->Initialize();
    S->Initialize();

    // take a single timestep
    double dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    pk->AdvanceStep("", "next");

    CHECK_CLOSE(1.0,
                (*S->Get<CompositeVector>("primary", "")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
    CHECK_CLOSE(2.0,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);

    pk->CommitStep("", "next");
    S->advance_time("", dt);
    S->advance_cycle("", 1);

    dt = pk->get_dt();
    S->advance_time("next", dt);
    S->advance_cycle("next", 1);
    pk->AdvanceStep("", "next");

    CHECK_CLOSE(2.0,
                (*S->Get<CompositeVector>("primary", "")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
    CHECK_CLOSE(3.0,
                (*S->Get<CompositeVector>("primary", "next")
                    .ViewComponent("cell", false))[0][0],
                1.e-10);
  }
}
