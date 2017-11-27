/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

/* Test basic implicit and explicit PKs.

At this point PKs manage memory and interface time integrators with the DAG.
These tests that functionality with a series of ODEs.
   
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

#include "PK_MixinMPC.hh"
#include "PK_MixinMPCAdvanceStepWeak.hh"
#include "PK_MixinMPCGetDtMin.hh"

#include "pks_test.hh"

using namespace Amanzi;

struct Run {
  Run(const Teuchos::RCP<State>& S_,
      const Teuchos::RCP<PK>& pk_,
      const Teuchos::RCP<TreeVectorSpace>& soln_map_)
      : S(S_),
        pk(pk_),
        soln_map(soln_map_) {}

  Teuchos::RCP<State> S;
  Teuchos::RCP<PK> pk;
  Teuchos::RCP<TreeVectorSpace> soln_map;
};

//
// Helper function that mocks a coordinator to run the test.
// ================================================================================
std::pair<double,double> run_test(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<PK>& pk,
        const Teuchos::RCP<TreeVectorSpace>& soln_map) {

  auto soln_old = Teuchos::rcp(new TreeVector(*soln_map, INIT_MODE_NOALLOC));
  auto soln_next = Teuchos::rcp(new TreeVector(*soln_map, INIT_MODE_NOALLOC));
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


//
// Creates an explicitly integrable PK
// ============================================================================
template<class MPC_t, class PK_A_t, class PK_B_t>
std::unique_ptr<Run>
createRun(const std::string& mpc_name, const std::string& pk_A_name, const std::string& pk_B_name) {
  std::cout << "Test: " << mpc_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

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


  // create the PKs/MPCs
  auto soln = Teuchos::rcp(new TreeVectorSpace());
  auto sv1 = Teuchos::rcp(new TreeVectorSpace());
  auto sv2 = Teuchos::rcp(new TreeVectorSpace());
  soln->PushBack(sv1);
  soln->PushBack(sv2);

  Teuchos::RCP<PK> pk_mpc;

  
  auto pk_tree_mpc = Teuchos::rcp(new Teuchos::ParameterList(mpc_name));
  auto mpc = Teuchos::rcp(new MPC_t(pk_tree_mpc, global_list, S, soln));

  auto pk_tree_A = Teuchos::rcp(new Teuchos::ParameterList(pk_A_name));
  auto pk_a = Teuchos::rcp(new PK_A_t(pk_tree_A, global_list, S, sv1));

  auto pk_tree_B = Teuchos::rcp(new Teuchos::ParameterList(pk_B_name));
  auto pk_b = Teuchos::rcp(new PK_B_t(pk_tree_B, global_list, S, sv2));
  auto sub_pks = std::vector<Teuchos::RCP<PK> >{pk_a, pk_b};
  mpc->SetChildren(sub_pks);

  return std::make_unique<Run>(S,mpc,soln);
}


SUITE(PKS_MPC) {

  // weak MPC coupling two FE PKs
  TEST(WEAK_BC_FORWARD_EULER) {
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> > PK_B_t;
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> > PK_C_t;
    typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default,PK> > > > MPC_t;

    auto run = createRun<MPC_t, PK_B_t, PK_C_t>("BC weak forward euler", "B, forward euler", "C, forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln_map);

    // check B soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.59374, (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 1.e-4);

    // check C soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 0.4);
    CHECK_CLOSE(2.33463, (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    
    CHECK_EQUAL(10, nsteps.first);
  }

  // weak MPC coupling one FE and one RK4 PKs
  TEST(WEAK_BC_FE_RK) {
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> > PK_B_t;
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> > PK_C_t;
    typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default,PK> > > > MPC_t;

    auto run = createRun<MPC_t, PK_B_t, PK_C_t>("BC weak mixed explicit", "B, RK4", "C, forward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln_map);


    // check B soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 1.e-5);
    // check C soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 0.4);
    CHECK_CLOSE(2.33463, (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    
    CHECK_EQUAL(10, nsteps.first);

  }

  // weak MPC coupling one FE and one implicit BDF with a fixed timestep 
  TEST(WEAK_BC_FE_BDF) {
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> > PK_B_t;
    typedef PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> > PK_C_t;
    typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default,PK> > > > MPC_t;

    auto run = createRun<MPC_t, PK_B_t, PK_C_t>("BC weak imex", "B, forward euler", "C, backward euler");
    auto nsteps = run_test(run->S, run->pk, run->soln_map);


    // check B soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.59374, (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 1.e-4);

    // check C soln
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 0.6);
    CHECK_CLOSE(3.27476584420779, (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 1.e-8);
    
    CHECK_EQUAL(10, nsteps.first);

  }


  // weak MPC coupling one FE and one implicit BDF with a variable timestep in
  // which the BDF DOES fail, forcing the explicit to back up
  TEST(WEAK_BC_FE_BDF_FAILING) {
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> > PK_B_t;
    typedef PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> > PK_C_t;
    typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default,PK> > > > MPC_t;

    auto run = createRun<MPC_t, PK_B_t, PK_C_t>("BC weak imex variable dt", "B, RK4, large step", "C, backward euler, large step");
    auto nsteps = run_test(run->S, run->pk, run->soln_map);


    // check B soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.71826, (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 1.e-4);

    // check C soln
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 0.6);
    CHECK_CLOSE(3.01043, (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    
    CHECK_EQUAL(189, nsteps.first);
    CHECK(nsteps.second > 0);

  }


  // weak MPC coupling one FE and one implicit BDF with a variable timestep in
  // which the BDF DOES fail, and subcycles to keep up with the explicit
  TEST(WEAK_BC_FE_BDF_FAILING2) {
    typedef PK_Adaptor<PK_ODE_Explicit<PK_MixinExplicit<PK_MixinLeaf<PK_Default> >, DudtEvaluatorB> > PK_B_t;
    typedef PK_Adaptor<PK_ODE_Implicit<PK_MixinImplicitSubcycled<PK_MixinLeaf<PK_Default> >, DudtEvaluatorC> > PK_C_t;
    typedef PK_Adaptor<PK_MixinMPCAdvanceStepWeak<PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default,PK> > > > MPC_t;

    auto run = createRun<MPC_t, PK_B_t, PK_C_t>("BC weak imex variable dt", "B, RK4", "C, backward euler, large step");
    auto nsteps = run_test(run->S, run->pk, run->soln_map);


    // check B soln
    CHECK_CLOSE(std::exp(1), (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 0.15);
    CHECK_CLOSE(2.71826, (*run->S->Get<CompositeVector>("primaryB")
                      .ViewComponent("cell",false))[0][0], 1.e-4);

    // check C soln
    CHECK_CLOSE(std::exp(1.0), (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 0.1);
    CHECK_CLOSE(2.77256, (*run->S->Get<CompositeVector>("primaryC")
                      .ViewComponent("cell",false))[0][0], 1.e-4);
    
    CHECK_EQUAL(10, nsteps.first);
  }

}
