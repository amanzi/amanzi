/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "EvaluatorPrimary.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Tag.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/* ******************************************************************
* Equation A = 2*B
****************************************************************** */
class AEvaluator
    : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
public:
  AEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist) {
    dependencies_.insert(std::make_pair(Key("fb"), Tags::DEFAULT));
    comp_ = plist.get<std::string>("component", "cell");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this));
  };

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override {
    auto& result_c = *results[0]->ViewComponent(comp_);
    const auto& fb_c = *S.Get<CompositeVector>("fb").ViewComponent(comp_);

    for (int c = 0; c != result_c.MyLength(); ++c) {
      result_c[0][c] = 2 * fb_c[0][c];
    }
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*> &results) override {
    auto& result_c = *results[0]->ViewComponent(comp_);

    if (wrt_key == "fb") {
      for (int c = 0; c != result_c.MyLength(); ++c) {
        result_c[0][c] = 2.0;
      }
    }
  }

  std::string comp_;
};


SUITE(EVALUATORS_CV) {
  TEST(PRIMARY_CV) {
    // Tests a primary variable evaluator
    std::cout << "\nPrimary Variable Test" << std::endl;
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName("fa");

    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector,CompositeVectorSpace>(
        "fa", Tags::DEFAULT, "fa", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    // Setup fields and marked as initialized
    S.Setup();
    S.GetW<CompositeVector>("fa", "fa").PutScalar(1.0);
    S.GetRecordW("fa", "fa").set_initialized();
    S.Initialize();

    // provides
    auto tag1 = make_tag("old");
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", Tags::DEFAULT));      // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("other", Tags::DEFAULT));  // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", tag1));              // self

    // dependencies -- none
    auto tag2 = make_tag("other");
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", Tags::DEFAULT));     // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", Tags::DEFAULT));  // not other
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", tag2));              // not other

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", Tags::DEFAULT));

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", Tags::DEFAULT));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fa", Tags::DEFAULT));

    // mark as changed
    auto eval = S.GetEvaluatorPtr("fa", Tags::DEFAULT);
    auto eval_p = Teuchos::rcp_dynamic_cast<
        EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", Tags::DEFAULT));

    CHECK_CLOSE(1.0, (*S.GetDerivative<CompositeVector>(
        "fa", Tags::DEFAULT, "fa", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);
  }


  TEST(SECONDARY_CV_DEFAULT_USAGE) {
    // Tests a secondary variable evaluator in standard usage
    std::cout << "\nSecondary Variable Test (default)" << std::endl;
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb")
        .SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
        
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)->SetGhosted(true)
                      ->SetComponent("cell", AmanziMesh::CELL, 1);

    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
        "fa", Tags::DEFAULT, "fb", Tags::DEFAULT);

    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
        "fa", Tags::DEFAULT, "fa", Tags::DEFAULT);

    // setup 
    S.Setup();

    // Check that info from fa made it into fb under the default data policy
    CHECK(S.Get<CompositeVector>("fa", Tags::DEFAULT).HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fa", Tags::DEFAULT).Mesh() == S.GetMesh("domain"));

    // initialize
    S.GetW<CompositeVector>("fb", Tags::DEFAULT, "fb").PutScalar(3.0);
    S.GetRecordW("fb", "fb").set_initialized();
    S.Initialize();

    // provides
    auto tag1 = make_tag("old");
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", Tags::DEFAULT));   // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", tag1));           // other tag
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fb", Tags::DEFAULT));  // other key

    // dependencies -- fb
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", Tags::DEFAULT));     // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", tag1));              // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", Tags::DEFAULT));  // not other
    CHECK(S.GetEvaluator("fa").IsDependency(S, "fb", Tags::DEFAULT));      // but fb is!

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", Tags::DEFAULT));

    // check the value and derivative
    CHECK_CLOSE(6.0, (*S.Get<CompositeVector>("fa", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);
    CHECK_CLOSE(2.0, (*S.GetDerivative<CompositeVector>(
        "fa", Tags::DEFAULT, "fb", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", Tags::DEFAULT));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fb", Tags::DEFAULT));

    // check the value and derivative are still the same
    CHECK_CLOSE(6.0, (*S.Get<CompositeVector>("fa", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);
    CHECK_CLOSE(2.0, (*S.GetDerivative<CompositeVector>(
        "fa", Tags::DEFAULT, "fb", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);

    // change the primary and mark as changed
    S.GetW<CompositeVector>("fb", Tags::DEFAULT, "fb").PutScalar(14.0);
    auto eval = S.GetEvaluatorPtr("fb");
    auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector,CompositeVectorSpace>>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    // now the value is different, and so it has changed
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", Tags::DEFAULT));

    // check the values
    CHECK_CLOSE(28.0, (*S.Get<CompositeVector>("fa", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);
    CHECK_CLOSE(2.0, (*S.GetDerivative<CompositeVector>(
        "fa", Tags::DEFAULT, "fb", Tags::DEFAULT).ViewComponent("cell"))[0][0], 1.0e-10);

    // check self-derivative is not working
    // CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", Tags::DEFAULT));
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES) {
    // Tests two secondaries depending upon one primary, making sure they correctly set structure
    std::cout << "\nSecondary Variable Test (joint deps)" << std::endl;
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb");
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)->SetGhosted(true)
                      ->SetComponent("cell", AmanziMesh::CELL, 1);

    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    // make another secondary
    Teuchos::ParameterList ec_list;
    ec_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ec_list.setName("fc");
    ec_list.set("tag", "");
    ec_list.set("component", "face");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fc")
        .SetMesh(mesh)->SetGhosted(true)
                      ->SetComponent("face", AmanziMesh::FACE, 2);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fc_eval = Teuchos::rcp(new AEvaluator(ec_list));
    S.SetEvaluator("fc", Tags::DEFAULT, fc_eval);
    
    // setup 
    S.Setup();

    CHECK(S.Get<CompositeVector>("fa", Tags::DEFAULT).Mesh() == S.GetMesh("domain"));

    // b gets both
    CHECK(S.Get<CompositeVector>("fb", Tags::DEFAULT).HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fb", Tags::DEFAULT).HasComponent("face"));

    // a still has just cell
    CHECK(S.Get<CompositeVector>("fa", Tags::DEFAULT).HasComponent("cell"));
    CHECK(!S.Get<CompositeVector>("fa", Tags::DEFAULT).HasComponent("face"));
    CHECK(S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT).HasComponent("cell"));
    CHECK(!S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT).HasComponent("face"));

    // c still has just face
    CHECK(!S.Get<CompositeVector>("fc", Tags::DEFAULT).HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fc", Tags::DEFAULT).HasComponent("face"));
    CHECK(!S.GetDerivative<CompositeVector>("fc", Tags::DEFAULT, "fb", Tags::DEFAULT).HasComponent("cell"));
    CHECK(S.GetDerivative<CompositeVector>("fc", Tags::DEFAULT, "fb", Tags::DEFAULT).HasComponent("face"));
  }

  
  TEST(SECONDARY_CV_THROWS_PARENT_CHILD_INCONSISTENT_MESH) {
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb").SetMesh(mesh2);
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    CHECK_THROW(S.Setup(), Errors::Message);
  }


  TEST(SECONDARY_CV_THROWS_PARENT_CHILD_INCONSISTENT_STRUCTURE) {
    // Test that Setup() throws if parent and child structure/components are different
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb")
        .SetMesh(mesh)
        ->SetComponent("face", AmanziMesh::FACE, 1);
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    CHECK_THROW(S.Setup(), Errors::Message);
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_INCONSISTENT_MESH) {
    // Test that Setup() throws if two parents have inconsistent meshes and try to both set their child
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb");
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);


    // make another secondary
    Teuchos::ParameterList ec_list;
    ec_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ec_list.setName("fc");
    ec_list.set("tag", "");
    ec_list.set("component", "face");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fc")
        .SetMesh(mesh2)
        ->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 2);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fc_eval = Teuchos::rcp(new AEvaluator(ec_list));
    S.SetEvaluator("fc", Tags::DEFAULT, fc_eval);
    
    // setup 
    CHECK_THROW(S.Setup(), Errors::Message);
  }



  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_INCONSISTENT_STRUCTURE) {
    // Test that Setup() throws if two parents have inconsistent structure and both try to set the parent
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb");
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 2);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);


    // make another secondary
    Teuchos::ParameterList ec_list;
    ec_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ec_list.setName("fc");
    ec_list.set("tag", "");
    ec_list.set("component", "face");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fc")
        .SetMesh(mesh2)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fc_eval = Teuchos::rcp(new AEvaluator(ec_list));
    S.SetEvaluator("fc", Tags::DEFAULT, fc_eval);
    
    // setup 
    CHECK_THROW(S.Setup(), Errors::Message);
  }
  

  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_REVERSED_DEPENDENCY) {
    // Test the reversed structure case, where a parent takes its structure from a child
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb")
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    ea_list.set("consistency policy", "take from child");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa");
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);
    
    // setup 
    S.Setup();

    CHECK(S.Get<CompositeVector>("fa", Tags::DEFAULT).Mesh() == S.GetMesh("domain"));
    CHECK(S.Get<CompositeVector>("fb", Tags::DEFAULT).HasComponent("cell"));
    CHECK(S.Get<CompositeVector>("fa", Tags::DEFAULT).HasComponent("cell"));
    CHECK(S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT).HasComponent("cell"));
  }


  TEST(SECONDARY_CV_JOINT_DEPENDENCIES_REVERSED_DEPENDENCY_WITHOUT_FLAG_THROWS) {
    // Test the reversed structure case, where a parent takes its structure from a child
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
    S.Require<CompositeVector,CompositeVectorSpace>("fb", Tags::DEFAULT, "fb")
        .SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(es_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // make the secondary
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
        .SetMesh(mesh);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);
    
    // setup 
    CHECK_THROW(S.Setup(), Errors::Message);
  }
}
