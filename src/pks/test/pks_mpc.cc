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

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
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

SUITE(PKS_MPC)
{
  // weak MPC coupling two FE PKs
  TEST(SEQUENTIAL_BC_FORWARD_EULER)
  {
    typedef PK_Explicit_Adaptor<
      PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorB>>
      PK_B_t;
    typedef PK_Explicit_Adaptor<
      PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
                      DudtEvaluatorC>>
      PK_C_t;
    typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
      PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>
      MPC_t;

    auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t>(
      "BC weak forward euler", "B, forward euler", "C, forward euler");
    auto nsteps = run_test(run->S, run->pk);

    // check B soln -- same as B_FORWARD_EULER
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.15);
    CHECK_CLOSE(2.59374, // calculated via test.py
                run->S->Get<CompositeVector>("primaryB")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);

    // check C soln -- same as C_FORWARD_EULER
    CHECK_CLOSE(std::exp(1),
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                0.4);
    CHECK_CLOSE(2.33463, // calculated via test.py
                run->S->Get<CompositeVector>("primaryC")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-4);

    CHECK_EQUAL(10, nsteps.first);
  }

  // // weak MPC coupling one FE and one RK4 PKs
  // TEST(SEQUENTIAL_BC_FE_RK) {
  //   typedef PK_Explicit_Adaptor<PK_ODE_Explicit<
  //       PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorB>>
  //       PK_B_t;
  //   typedef PK_Explicit_Adaptor<PK_ODE_Explicit<
  //       PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorC>>
  //       PK_C_t;
  //   typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
  //       PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>
  //       MPC_t;

  //   auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t>(
  //       "BC weak mixed explicit", "B, RK4", "C, forward euler");
  //   auto nsteps = run_test(run->S, run->pk);

  //   // check B soln - same as B_RK4
  //   CHECK_CLOSE(std::exp(1),
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-5);
  //   // check C soln - same as C_FORWARD_EULER
  //   CHECK_CLOSE(std::exp(1),
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.4);
  //   CHECK_CLOSE(2.33463,
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);

  //   CHECK_EQUAL(10, nsteps.first);
  //   CHECK_EQUAL(0, nsteps.second);
  // }

  // // weak MPC coupling one FE and one implicit BDF with a fixed timestep
  // TEST(SEQUENTIAL_BC_FE_BDF) {
  //   typedef PK_Explicit_Adaptor<PK_ODE_Explicit<
  //       PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorB>>
  //       PK_B_t;
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorC>>
  //       PK_C_t;
  //   typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
  //       PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>
  //       MPC_t;

  //   auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t>(
  //       "BC weak imex", "B, forward euler", "C, backward euler");
  //   auto nsteps = run_test(run->S, run->pk);

  //   // check B soln - same as B_FORWARD_EULER
  //   CHECK_CLOSE(std::exp(1),
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.15);
  //   CHECK_CLOSE(2.59374,
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);

  //   // check C soln - same as C_BACKWARD_EULER
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.6);
  //   CHECK_CLOSE(3.27476584420779,
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-8);

  //   CHECK_EQUAL(10, nsteps.first);
  //   CHECK_EQUAL(0, nsteps.second);
  // }

  // // weak MPC coupling one FE and one implicit BDF with a variable timestep
  // in
  // // which the BDF DOES fail, forcing the explicit to back up
  // TEST(SEQUENTIAL_BC_FE_BDF_FAILING) {
  //   typedef PK_Explicit_Adaptor<PK_ODE_Explicit<
  //       PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorB>>
  //       PK_B_t;
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorC>>
  //       PK_C_t;
  //   typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
  //       PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>
  //       MPC_t;

  //   auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t>(
  //       "BC weak imex variable dt", "B, RK4, large step",
  //       "C, backward euler, large step");
  //   auto nsteps = run_test(run->S, run->pk);

  //   // check B soln - no analogue
  //   CHECK_CLOSE(std::exp(1),
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.15);
  //   CHECK_CLOSE(2.71826,
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);

  //   // check C soln
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.6);
  //   CHECK_CLOSE(3.02734,
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);
  //   CHECK_EQUAL(69, nsteps.first);
  //   CHECK_EQUAL(2, nsteps.second);
  // }

  // // weak MPC coupling one FE and one implicit BDF with a variable timestep
  // in
  // // which the BDF DOES fail, and subcycles to keep up with the explicit
  // TEST(SEQUENTIAL_BC_FE_BDF_FAILING2) {
  //   typedef PK_Explicit_Adaptor<PK_ODE_Explicit<
  //       PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorB>>
  //       PK_B_t;
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorC>>
  //       PK_C_t;
  //   typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
  //       PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>
  //       MPC_t;

  //   auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t>(
  //       "BC weak imex variable dt", "B, RK4", "C, backward euler, large
  //       step");
  //   auto nsteps = run_test(run->S, run->pk);

  //   // check B soln - note this is the same as pks_ode:B_RK4
  //   CHECK_CLOSE(std::exp(1),
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-5);

  //   // check C soln - note this is the same as
  //   // pks_ode:C_BACKWARD_EULER_SUBCYCLED_MULTIPLE
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.1);
  //   CHECK_CLOSE(2.79649,
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);

  //   CHECK_EQUAL(10, nsteps.first);
  //   CHECK_EQUAL(0, nsteps.second);
  // }

  // // Globally implicit, fixed timestep
  // TEST(IMPLICIT_BC) {
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorB>>
  //       PK_B_t;
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorC>>
  //       PK_C_t;
  //   typedef PK_Implicit_Adaptor<PK_MixinImplicit<
  //       PK_MixinMPCImplicit<PK_Default, PK_Implicit<TreeVector>>>>
  //       MPC_t;
  //   auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t, PK_Implicit<>>(
  //       "BC global implicit", "B, backward euler", "C, backward euler");
  //   auto nsteps = run_test(run->S, run->pk);

  //   // check B soln, same as B_BACKWARD_EULER
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.15);
  //   CHECK_CLOSE(2.867971990790009,
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-8);

  //   // check C soln, same as C_BACKWARD_EULER
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.6);
  //   CHECK_CLOSE(3.27476584420779,
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-8);

  //   CHECK_EQUAL(10, nsteps.first);
  //   CHECK_EQUAL(0, nsteps.second);
  // }

  // // Globally implicit, variable timestep
  // TEST(IMPLICIT_BC_VARIABLE_TS) {
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorB>>
  //       PK_B_t;
  //   typedef PK_Implicit_Adaptor<PK_ODE_Implicit<
  //       PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>,
  //       DudtEvaluatorC>>
  //       PK_C_t;
  //   typedef PK_Implicit_Adaptor<PK_MixinImplicit<
  //       PK_MixinMPCImplicit<PK_Default, PK_Implicit<TreeVector>>>>
  //       MPC_t;
  //   auto run = createRunMPC<MPC_t, PK_B_t, PK_C_t, PK_Implicit<>>(
  //       "BC global implicit variable", "B, backward euler",
  //       "C, backward euler");
  //   auto nsteps = run_test(run->S, run->pk);

  //   // no analogue for either of these
  //   // check B soln
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.04);
  //   CHECK_CLOSE(2.74863,
  //               run->S->Get<CompositeVector>("primaryB")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);

  //   // check C soln
  //   CHECK_CLOSE(std::exp(1.0),
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               0.08);
  //   CHECK_CLOSE(2.79308,
  //               run->S->Get<CompositeVector>("primaryC")
  //                     .ViewComponent<MirrorHost>("cell", false)(0,0),
  //               1.e-4);

  //   CHECK_EQUAL(96, nsteps.first);
  // }
}
