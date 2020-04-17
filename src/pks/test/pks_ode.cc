/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/* Test basic implicit and explicit PKs.

At this point PKs manage memory and interface time integrators with the DAG.
These tests that functionality with a series of ODEs.

*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

#include "PK.hh"
#include "PK_Adaptors.hh"
#include "PK_Default.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinExplicitSubcycled.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinImplicitSubcycled.hh"
#include "PK_MixinLeaf.hh"
#include "PK_MixinPredictorCorrector.hh"

#include "pks_test_harness.hh"
#include "test_pks.hh"

using namespace Amanzi;

//
// Creates an explicitly integrable PK
// ============================================================================
std::unique_ptr<Run>
createExplicit(const std::string& eqn_name, const std::string& ti_name)
{
  std::string pk_name = eqn_name + ", " + ti_name;
  std::cout << std::endl
            << "Test: " << pk_name << std::endl
            << "================================================================================" << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<PK_Explicit<>> pk;
  if (eqn_name == "A") {
    pk =
      Teuchos::rcp(new PK_Explicit_Adaptor<PK_ODE_Explicit<
                     PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                     DudtEvaluatorA>>(pk_tree, global_list, S));
  } else if (eqn_name == "B") {
    pk =
      Teuchos::rcp(new PK_Explicit_Adaptor<PK_ODE_Explicit<
                     PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                     DudtEvaluatorB>>(pk_tree, global_list, S));
  } else if (eqn_name == "C") {
    pk =
      Teuchos::rcp(new PK_Explicit_Adaptor<PK_ODE_Explicit<
                     PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                     DudtEvaluatorC>>(pk_tree, global_list, S));
  } else {
    AMANZI_ASSERT(false);
  }

  return std::make_unique<Run>(S, pk);
}

//
// Creates an explicit subcycled PK
// ============================================================================
std::unique_ptr<Run>
createExplicitSubcycled(const std::string& eqn_name, const std::string& ti_name)
{
  std::string pk_name = eqn_name + ", " + ti_name + " subcycled";
  std::cout << std::endl
            << "Test: " << pk_name << std::endl
            << "================================================================================" << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<PK_Explicit<>> pk;
  if (eqn_name == "A") {
    pk = Teuchos::rcp(
      new PK_Explicit_Adaptor<PK_ODE_Explicit<
        PK_MixinExplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorA>>(pk_tree, global_list, S));
  } else if (eqn_name == "B") {
    pk = Teuchos::rcp(
      new PK_Explicit_Adaptor<PK_ODE_Explicit<
        PK_MixinExplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorB>>(pk_tree, global_list, S));
  } else if (eqn_name == "C") {
    pk = Teuchos::rcp(
      new PK_Explicit_Adaptor<PK_ODE_Explicit<
        PK_MixinExplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorC>>(pk_tree, global_list, S));
  } else {
    AMANZI_ASSERT(false);
  }

  return std::make_unique<Run>(S, pk);
}

//
// Creates an implicitly integrable PK
// ============================================================================
std::unique_ptr<Run>
createImplicit(const std::string& eqn_name, const std::string& qualifier = "")
{
  std::string pk_name = eqn_name + ", backward euler";
  if (!qualifier.empty()) pk_name = pk_name + ", " + qualifier;
  std::cout << std::endl
            << "Test: " << pk_name << std::endl
            << "================================================================================" << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<PK_Implicit<>> pk;
  if (eqn_name == "A") {
    pk =
      Teuchos::rcp(new PK_Implicit_Adaptor<PK_ODE_Implicit<
                     PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                     DudtEvaluatorA>>(pk_tree, global_list, S));
  } else if (eqn_name == "B") {
    pk =
      Teuchos::rcp(new PK_Implicit_Adaptor<PK_ODE_Implicit<
                     PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                     DudtEvaluatorB>>(pk_tree, global_list, S));
  } else if (eqn_name == "C") {
    pk =
      Teuchos::rcp(new PK_Implicit_Adaptor<PK_ODE_Implicit<
                     PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                     DudtEvaluatorC>>(pk_tree, global_list, S));
  } else {
    AMANZI_ASSERT(false);
  }

  return std::make_unique<Run>(S, pk);
}

// Creates a subcycled implicit PK.  These don't fail, they simply subcycle as
// needed to match the full timestep.
// ============================================================================
std::unique_ptr<Run>
createImplicitSubcycled(const std::string& eqn_name,
                        const std::string& qualifier = "")
{
  std::string pk_name = eqn_name + ", backward euler subcycled";
  if (!qualifier.empty()) pk_name = pk_name + ", " + qualifier;

  std::cout << std::endl
            << "Test: " << pk_name << std::endl
            << "================================================================================" << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));
  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<PK_Implicit<>> pk;
  if (eqn_name == "A") {
    pk = Teuchos::rcp(
      new PK_Implicit_Adaptor<PK_ODE_Implicit<
        PK_MixinImplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorA>>(pk_tree, global_list, S));
  } else if (eqn_name == "B") {
    pk = Teuchos::rcp(
      new PK_Implicit_Adaptor<PK_ODE_Implicit<
        PK_MixinImplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorB>>(pk_tree, global_list, S));
  } else if (eqn_name == "C") {
    pk = Teuchos::rcp(
      new PK_Implicit_Adaptor<PK_ODE_Implicit<
        PK_MixinImplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorC>>(pk_tree, global_list, S));
  } else {
    AMANZI_ASSERT(false);
  }

  return std::make_unique<Run>(S, pk);
}

SUITE(PKS_ODE)
{
  // Forward Euler tests with each of 3 PKs
  TEST(A_FORWARD_EULER)
  {
    auto run = createExplicit("A", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2.0,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-10);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_FORWARD_EULER)
  {
    auto run = createExplicit("B", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.15);
    CHECK_CLOSE(2.59374, // calculated via test.py
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_FORWARD_EULER)
  {
    auto run = createExplicit("C", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.4);
    CHECK_CLOSE(2.33463, // calculated via test.py
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  // Runge Kutta (multistage explicit) tests with each of 3 PKs
  TEST(A_RK4)
  {
    auto run = createExplicit("A", "RK4");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_RK4)
  {
    auto run = createExplicit("B", "RK4");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-5);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_RK4)
  {
    auto run = createExplicit("C", "RK4");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  // Subcycled Forward Euler tests with each of 3 PKs.  Note that this answer
  // is identical to the Forward Euler -- it takes 10 subcycled steps for an
  // outer step that is 10 times as big.
  TEST(A_FORWARD_EULER_SUBCYCLED)
  {
    auto run = createExplicitSubcycled("A", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2.0,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-10);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_FORWARD_EULER_SUBCYCLED)
  {
    auto run = createExplicitSubcycled("B", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.15);
    CHECK_CLOSE(2.59374,
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_FORWARD_EULER_SUBCYCLED)
  {
    auto run = createExplicitSubcycled("C", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.4);
    CHECK_CLOSE(2.33463,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(A_RK4_SUBCYCLED)
  {
    auto run = createExplicitSubcycled("A", "RK4");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2.0,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-10);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_RK4_SUBCYCLED)
  {
    auto run = createExplicitSubcycled("B", "RK4");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-5);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_RK4_SUBCYCLED)
  {
    auto run = createExplicitSubcycled("C", "RK4");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  // Implicit (single stage)
  TEST(A_BACKWARD_EULER)
  {
    auto run = createImplicit("A");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_BACKWARD_EULER)
  {
    auto run = createImplicit("B");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.15);
    CHECK_CLOSE(2.867971990790009,
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER)
  {
    auto run = createImplicit("C");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.6);
    CHECK_CLOSE(3.27476584420779,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER_VARIABLE)
  {
    auto run = createImplicit("C", "large step");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.6);
    CHECK_CLOSE(3.02734,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(69, nsteps.first);
    CHECK_EQUAL(2, nsteps.second);
  }

  // Implicit (subcycled)
  TEST(A_BACKWARD_EULER_SUBCYCLED)
  {
    auto run = createImplicitSubcycled("A");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-8);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(B_BACKWARD_EULER_SUBCYCLED)
  {
    auto run = createImplicitSubcycled("B");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.15);
    CHECK_CLOSE(2.86608,
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  TEST(C_BACKWARD_EULER_SUBCYCLED)
  {
    // note this is identical to C_BACKWARD_EULER_VARIABLE
    auto run = createImplicitSubcycled("C");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.6);
    CHECK_CLOSE(3.02734,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
    CHECK_EQUAL(69,
                run->S->Get<int>(
                  "cycle", "C, backward euler subcycled implicit inner next"));
  }

  TEST(C_BACKWARD_EULER_SUBCYCLED_MULTIPLE)
  {
    // force a subcyle
    auto run = createImplicitSubcycled("C", "forced");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.6);
    CHECK_CLOSE(2.79649,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);

    // this is not the total count, but the count of the last outer step's inner
    // steps
    CHECK_EQUAL(
      17,
      run->S->Get<int>(
        "cycle", "C, backward euler subcycled, forced implicit inner next"));
  }

  TEST(C_PREDICTOR_CORRECTOR)
  {
    using PK_t = PK_ImplicitExplicit_Adaptor<PK_ODE_Implicit<
      PK_ODE_Explicit<
        PK_MixinPredictorCorrector<PK_MixinLeafCompositeVector<PK_Default>>,
        DudtEvaluatorC>,
      DudtEvaluatorC>>;
    auto run = createRunODE<PK_t>("C predictor corrector");
    auto nsteps = run_test(run->S, run->pk);

    // note, worse error, but fewer timesteps, than the implicit version with
    // linear extrapolation
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.6);
    CHECK_CLOSE(3.02976,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);
    CHECK_EQUAL(65, nsteps.first);
    CHECK_EQUAL(2, nsteps.second);
  }
}
