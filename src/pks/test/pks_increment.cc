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

#include "AmanziTypes.hh"
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
#include "test_increment_pk.hh"
#include "test_pks.hh"

using namespace Amanzi;

std::unique_ptr<Run>
create(const std::string& eqn_name)
{
  std::string pk_name = eqn_name;
  std::cout << "Test: " << pk_name << std::endl;

  auto global_list = Teuchos::getParametersFromXmlFile("test/pks_ode.xml");
  auto pk_tree = Teuchos::rcp(new Teuchos::ParameterList(pk_name));

  auto S = Teuchos::rcp(new State(global_list->sublist("state")));

  auto comm = getDefaultComm();

  // create mesh
  Teuchos::ParameterList& regions_list = global_list->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
    meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<PK> pk = Teuchos::rcp(
    new PK_Adaptor<PK_Veg<PK_MixinLeaf<PK_Default, PFTList, PFTListSpace>>>(
      pk_tree, global_list, S));
  return std::make_unique<Run>(S, pk);
}

SUITE(PKS_INCREMENT)
{
  //
  // Creates an increment PK, i.e. a PK that does its own advance and cannot
  // fail
  // ============================================================================
  TEST(PK_VEG)
  {
    auto run = create("Veg");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(0.1, run->S->Get<PFTList>("PFT_A")[0].Bleaf, 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }

  //
  // Create an incrememt PK, but couple it to a subcycled implicit PK.  The
  // subcycled implicit PK can and does fail, but the increment PK cannot fail
  // (it might have internal state, ugh).
  // =============================================================================
  TEST(PK_VEG_WATER)
  {
    using PK_Veg_t =
      PK_Adaptor<PK_Veg<PK_MixinLeaf<PK_Default, PFTList, PFTListSpace>>>;
    using PK_Water_t = PK_Implicit_Adaptor<PK_ODE_Implicit<
      PK_MixinImplicitSubcycled<PK_MixinLeafCompositeVector<PK_Default>>,
      DudtEvaluatorC>>;
    using MPC_t = PK_Adaptor<PK_MixinMPCAdvanceStepWeak<
      PK_MixinMPCGetDtMin<PK_MixinMPC<PK_Default, PK>>>>;

    auto run = createRunMPC<MPC_t, PK_Veg_t, PK_Water_t>(
      "BC weak mixed explicit", "Veg", "C, backward euler subcycled, forced");
    auto nsteps = run_test(run->S, run->pk);

    CHECK_CLOSE(0.1, run->S->Get<PFTList>("PFT_A")[0].Bleaf, 1.e-8);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
}
