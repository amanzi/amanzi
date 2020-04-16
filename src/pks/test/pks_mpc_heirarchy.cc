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
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Adaptors.hh"
#include "PK_Default.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinExplicitSubcycled.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinImplicitSubcycled.hh"
#include "PK_MixinLeaf.hh"

#include "PK_MixinMPC.hh"
#include "PK_MixinMPCAdvanceStepWeak.hh"
#include "PK_MixinMPCGetDtMin.hh"
#include "PK_MixinMPCImplicit.hh"

#include "pks_test_harness.hh"
#include "test_pks.hh"

using namespace Amanzi;

//
// Creates an MPC heirarchy
//
//              MPC_A
//             /     \
//         MPC_B      PK_A
//        /     \
//     PK_B    PK_C
// ============================================================================
template <class MPC_A_t, class MPC_A_Children_t, class PK_A_t, class MPC_B_t,
          class MPC_B_Children_t, class PK_B_t, class PK_C_t>
std::unique_ptr<Run>
createRun(const std::string& mpc_A_name, const std::string& pk_A_name,
          const std::string& mpc_B_name, const std::string& pk_B_name,
          const std::string& pk_C_name)
{
  std::cout << std::endl
            << "           Test:" << std::endl
            << "               " << mpc_A_name << std::endl
            << "             /                         \\" << std::endl
            << "   " << mpc_B_name << "            " << pk_A_name << std::endl
            << "   /                   \\" << std::endl
            << pk_C_name << "     " << pk_B_name << std::endl
            << "======================================================================" << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  // create the PKs/MPCs, in the same order that they will be constructed when
  // done recursively through the factory
  auto pk_tree_mpc_A = Teuchos::rcp(new Teuchos::ParameterList(mpc_A_name));
  auto mpc_a = Teuchos::rcp(new MPC_A_t(pk_tree_mpc_A, global_list, S));

  auto pk_tree_A = Teuchos::rcp(new Teuchos::ParameterList(pk_A_name));
  auto pk_a = Teuchos::rcp(new PK_A_t(pk_tree_A, global_list, S));

  auto pk_tree_mpc_B = Teuchos::rcp(new Teuchos::ParameterList(mpc_B_name));
  auto mpc_b = Teuchos::rcp(new MPC_B_t(pk_tree_mpc_B, global_list, S));
  mpc_a->SetChildren(
    std::vector<Teuchos::RCP<MPC_A_Children_t>>{ pk_a, mpc_b });

  auto pk_tree_B = Teuchos::rcp(new Teuchos::ParameterList(pk_B_name));
  auto pk_b = Teuchos::rcp(new PK_B_t(pk_tree_B, global_list, S));

  auto pk_tree_C = Teuchos::rcp(new Teuchos::ParameterList(pk_C_name));
  auto pk_c = Teuchos::rcp(new PK_C_t(pk_tree_C, global_list, S));
  mpc_b->SetChildren(std::vector<Teuchos::RCP<MPC_B_Children_t>>{ pk_b, pk_c });

  return std::make_unique<Run>(S, mpc_a);
}

SUITE(PKS_MPC_THREE)
{
  // Sequential coupling of all, min dt
  TEST(SEQUENTIAL_ABC)
  {
    using PK_A_t = PK_Explicit_Adaptor<
      PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorA>>;
    using PK_B_t = PK_Explicit_Adaptor<
      PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorB>>;
    using PK_C_t = PK_Explicit_Adaptor<
      PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorC>>;
    using MPC_t = PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
      PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>;

    auto run = createRun<MPC_t, PK, PK_A_t, MPC_t, PK, PK_B_t, PK_C_t>(
        "ABC weak forward euler",
      "A, forward euler",
      "BC weak forward euler",
      "B, forward euler",
      "C, forward euler");
    auto nsteps = run_test(run->S, run->pk);

    // check A soln -- same as A_FORWARD_EULER
    CHECK_CLOSE(2.0,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                1.e-10);

    // check B soln -- same as B_FORWARD_EULER
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                0.15);
    CHECK_CLOSE(2.59374,
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                1.e-4);

    // check C soln -- same as C_FORWARD_EULER
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                0.4);
    CHECK_CLOSE(2.33463,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                1.e-4);

    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  // sequentially couple A to implicitly coupled B&C.  Min dt
  TEST(SEQUENTIAL_A_IMPLICIT_BC)
  {
    using PK_A_t = PK_Explicit_Adaptor<
      PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorA>>;
    using PK_B_t = PK_Implicit_Adaptor<
      PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorB>>;
    using PK_C_t = PK_Implicit_Adaptor<
      PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorC>>;
    using MPC_A_t = PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
      PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>;
    using MPC_B_t = PK_Implicit_Adaptor<PK_MixinImplicit<
      PK_MixinMPCImplicit<PK_Default, PK_Implicit<TreeVector>>>>;

    auto run = createRun<MPC_A_t,
                         PK,
                         PK_A_t,
                         MPC_B_t,
                         PK_Implicit<TreeVector>,
                         PK_B_t,
                         PK_C_t>("ABC weak forward euler",
                                 "A, forward euler",
                                 "BC global implicit variable",
                                 "B, backward euler",
                                 "C, backward euler");
    auto nsteps = run_test(run->S, run->pk);

    // check A soln
    {
      CHECK_CLOSE(2.0,
                  run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                  1.e-10);
    }

    // check B soln -- same as IMPLICIT_BC_VARIABLE_TS
    {
      CHECK_CLOSE(std::exp(1.0),
                  run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                  0.04);
    }
    {
      CHECK_CLOSE(2.74863,
                  run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                  1.e-4);
    }

    // check C soln -- same as IMPLICIT_BC_VARIABLE_TS
    {
      CHECK_CLOSE(std::exp(1.0),
                  run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                  0.08);
    }
    {
      CHECK_CLOSE(2.79308,
                  run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                  1.e-4);
    }

    CHECK_EQUAL(96, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  // implicitly couple everything heirarchically
  TEST(IMPLICIT_ABC)
  {
    using PK_A_t = PK_Implicit_Adaptor<
      PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorA>>;
    using PK_B_t = PK_Implicit_Adaptor<
      PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorB>>;
    using PK_C_t = PK_Implicit_Adaptor<
      PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorC>>;
    using MPC_t = PK_Implicit_Adaptor<PK_MixinImplicit<
      PK_MixinMPCImplicit<PK_Default, PK_Implicit<TreeVector>>>>;

    auto run = createRun<MPC_t,
                         PK_Implicit<TreeVector>,
                         PK_A_t,
                         MPC_t,
                         PK_Implicit<TreeVector>,
                         PK_B_t,
                         PK_C_t>("BC global implicit variable",
                                 "A, backward euler",
                                 "BC global implicit variable",
                                 "B, backward euler",
                                 "C, backward euler");
    auto nsteps = run_test(run->S, run->pk);

    // check A soln
    CHECK_CLOSE(2.0,
                run->S->Get<CompositeVector>("primaryA")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                1.e-10);

    // check B soln -- same as IMPLICIT_BC_VARIABLE_TS
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                0.04);
    CHECK_CLOSE(2.74863,
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                1.e-4);

    // check C soln -- same as IMPLICIT_BC_VARIABLE_TS
    CHECK_CLOSE(std::exp(1.0),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                0.08);
    CHECK_CLOSE(2.79308,
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<DefaultHost>("cell", false)(0, 0),
                1.e-4);

    CHECK_EQUAL(96, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
}
