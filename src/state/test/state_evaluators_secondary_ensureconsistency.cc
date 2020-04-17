/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
  Extensive testing for secondary variable evaluators "EnsureConsistency"
  method, which propagates vector meta-data up and down the DAG.

  Authors: Ethan Coon
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorPrimary.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/* ******************************************************************
 * Equation A = 2*B
 ****************************************************************** */
class AEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  AEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.emplace_back(std::make_pair(Key("fb"), Key("")));
    comp_ = plist.get<std::string>("component", "cell");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AEvaluator(*this));
  };

  virtual void Evaluate_(const State& S,
                         const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->ViewComponent(comp_);
    const auto fb_c = S.Get<CompositeVector>("fb").ViewComponent(comp_);

    Kokkos::parallel_for(result_c.extent(0), KOKKOS_LAMBDA(const int c) {
      result_c(c, 0) = 2 * fb_c(c, 0);
    });
  }

  virtual void EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Key& wrt_tag,
    const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->ViewComponent(comp_);

    if (wrt_key == "fb") {
      Kokkos::parallel_for(result_c.extent(0),
                           KOKKOS_LAMBDA(const int c) { result_c(c, 0) = 2.; });
    }
  }

  std::string comp_;
};

SUITE(EVALUATORS_CV)
{
  TEST(PRIMARY_CV)
  {
    // Tests a primary variable evaluator
    auto comm = Amanzi::getDefaultComm();

    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fa");

    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fa", "");
    auto fa_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fa", fa_eval);

    // Setup fields and marked as initialized
    S.Setup();
    S.GetW<CompositeVector>("fa", "fa").putScalar(1.0);
    S.GetRecordW("fa", "fa").set_initialized();
    S.Initialize();

    // provides
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", ""));     // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("other", "")); // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", "old")); // self

    // dependencies -- none
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", ""));      // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", ""));   // not other
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", "other")); // not other

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fa", ""));

    // mark as changed
    auto eval = S.GetEvaluatorPtr("fa", "");
    auto eval_p = Teuchos::rcp_dynamic_cast<
      EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));
  }


  TEST(SECONDARY_CV_DEFAULT_USAGE)
  {
    // Tests a secondary variable evaluator in standard usage
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    // make the primary.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb");

    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);

    // setup
    S.Setup();

    // Check that info from fa made it into fb under the default data policy
    CHECK(S.Get<CompositeVector>("fa", "").HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fa", "").Mesh() == S.GetMesh("domain"));

    // initialize
    S.GetW<CompositeVector>("fb", "", "fb").putScalar(3.0);
    S.GetRecordW("fb", "fb").set_initialized();
    S.Initialize();

    // provides
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", ""));     // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", "old")); // other tag
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fb", ""));    // other key

    // dependencies -- fb
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", ""));    // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", "old")); // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", "")); // not other
    CHECK(S.GetEvaluator("fa").IsDependency(S, "fb", ""));     // but fb is!

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", ""));

    // check the value and derivative
    CHECK_CLOSE(
      6.0,
      S.Get<CompositeVector>("fa", "").ViewComponent<MirrorHost>(
        "cell", false)(0, 0),
      1.e-10);
    CHECK_CLOSE(2.0,
                S.GetDerivative<CompositeVector>("fa", "", "fb", "")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-10);

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", ""));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fb", ""));

    // check the value and derivative are still the same
    CHECK_CLOSE(
      6.0,
      S.Get<CompositeVector>("fa", "").ViewComponent<MirrorHost>(
        "cell", false)(0, 0),
      1.e-10);
    CHECK_CLOSE(2.0,
                S.GetDerivative<CompositeVector>("fa", "", "fb", "")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-10);

    // change the primary and mark as changed
    S.GetW<CompositeVector>("fb", "", "fb").putScalar(14.0);
    auto eval = S.GetEvaluatorPtr("fb");
    auto eval_p = Teuchos::rcp_dynamic_cast<
      EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    // now the value is different, and so it has changed
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", ""));

    // check the values
    CHECK_CLOSE(
      28.0,
      S.Get<CompositeVector>("fa", "").ViewComponent<MirrorHost>(
        "cell", false)(0, 0),
      1.e-10);
    CHECK_CLOSE(2.0,
                S.GetDerivative<CompositeVector>("fa", "", "fb", "")
                  .ViewComponent<MirrorHost>("cell", false)(0, 0),
                1.e-10);
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES)
  {
    // Tests two secondaries depending upon one primary, making sure they
    // correctly set structure
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb");
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);


    // make another secondary
    Teuchos::ParameterList ec_list;
    ec_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ec_list.setName("fc");
    ec_list.set("tag", "");
    ec_list.set("component", "face");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", "", "fc")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 2);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fc", "", "fb", "");
    auto fc_eval = Teuchos::rcp(new AEvaluator(ec_list));
    S.SetEvaluator("fc", fc_eval);

    // setup
    S.Setup();

    CHECK(S.Get<CompositeVector>("fa", "").Mesh() == S.GetMesh("domain"));

    // b gets both
    CHECK(S.Get<CompositeVector>("fb", "").HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fb", "").HasComponent("face"));

    // a still has just cell
    CHECK(S.Get<CompositeVector>("fa", "").HasComponent("cell"));
    CHECK(!S.Get<CompositeVector>("fa", "").HasComponent("face"));
    CHECK(S.GetDerivative<CompositeVector>("fa", "", "fb", "")
            .HasComponent("cell"));
    CHECK(!S.GetDerivative<CompositeVector>("fa", "", "fb", "")
             .HasComponent("face"));

    // c still has just face
    CHECK(!S.Get<CompositeVector>("fc", "").HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fc", "").HasComponent("face"));
    CHECK(!S.GetDerivative<CompositeVector>("fc", "", "fb", "")
             .HasComponent("cell"));
    CHECK(S.GetDerivative<CompositeVector>("fc", "", "fb", "")
            .HasComponent("face"));
  }


  TEST(SECONDARY_CV_THROWS_PARENT_CHILD_INCONSISTENT_MESH)
  {
    // Test that Setup() throws if parent and child meshes are different
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto mesh2 = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // check throws on mesh inconsistency
    State S;
    S.RegisterDomainMesh(mesh);
    S.RegisterMesh("m2", mesh2);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb")
      .SetMesh(mesh2);
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);

    CHECK_THROW(S.Setup(), Errors::Message);
  }

  TEST(SECONDARY_CV_THROWS_PARENT_CHILD_INCONSISTENT_STRUCTURE)
  {
    // Test that Setup() throws if parent and child structure/components are
    // different
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto mesh2 = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // check throws on mesh inconsistency
    State S;
    S.RegisterDomainMesh(mesh);
    S.RegisterMesh("m2", mesh2);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb")
      .SetMesh(mesh)
      ->SetComponent("face", AmanziMesh::FACE, 1);
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);

    CHECK_THROW(S.Setup(), Errors::Message);
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_INCONSISTENT_MESH)
  {
    // Test that Setup() throws if two parents have inconsistent meshes and try
    // to both set their child
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto mesh2 = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);
    S.RegisterMesh("m2", mesh2);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb");
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);


    // make another secondary
    Teuchos::ParameterList ec_list;
    ec_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ec_list.setName("fc");
    ec_list.set("tag", "");
    ec_list.set("component", "face");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", "", "fc")
      .SetMesh(mesh2)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 2);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fc", "", "fb", "");
    auto fc_eval = Teuchos::rcp(new AEvaluator(ec_list));
    S.SetEvaluator("fc", fc_eval);

    // setup
    CHECK_THROW(S.Setup(), Errors::Message);
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_INCONSISTENT_STRUCTURE)
  {
    // Test that Setup() throws if two parents have inconsistent structure and
    // both try to set the parent
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    auto mesh2 = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);
    S.RegisterMesh("m2", mesh2);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb");
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);


    // make another secondary
    Teuchos::ParameterList ec_list;
    ec_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ec_list.setName("fc");
    ec_list.set("tag", "");
    ec_list.set("component", "face");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", "", "fc")
      .SetMesh(mesh2)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fc", "", "fb", "");
    auto fc_eval = Teuchos::rcp(new AEvaluator(ec_list));
    S.SetEvaluator("fc", fc_eval);

    // setup
    CHECK_THROW(S.Setup(), Errors::Message);
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_REVERSED_DEPENDENCY)
  {
    // Test the reversed structure case, where a parent takes its structure from
    // a child
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    ea_list.set("consistency policy", "take from child");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa");
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);

    // setup
    S.Setup();

    CHECK(S.Get<CompositeVector>("fa", "").Mesh() == S.GetMesh("domain"));
    CHECK(S.Get<CompositeVector>("fb", "").HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fa", "").HasComponent("cell"));
    CHECK(S.GetDerivative<CompositeVector>("fa", "", "fb", "")
            .HasComponent("cell"));
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_REVERSED_DEPENDENCY_WITHOUT_FLAG_THROWS)
  {
    // Test the reversed structure case, where a parent takes its structure from
    // a child
    std::cout << "Secondary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    // make the primary
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<CompositeVector, CompositeVectorSpace>("fb", "", "fb")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fb_eval = Teuchos::rcp(
      new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);

    // setup
    CHECK_THROW(S.Setup(), Errors::Message);
  }
}
