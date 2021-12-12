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

#include "EvaluatorIndependent.hh"
#include "EvaluatorPrimary.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/* ******************************************************************
 * Equation A = 2*B
 ****************************************************************** */
class AEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  AEvaluator(Teuchos::ParameterList &plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.emplace_back(std::make_pair(Key{"fb"}, Key{""}));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this));
  };

  virtual void Evaluate_(const State &S, const std::vector<double*> &results) override {
    auto &fb = S.Get<double>("fb");
    (*results[0]) = 2 * fb;
  }

  virtual void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                          const Key &wrt_tag,
                                          const std::vector<double*> &results) override {
    if (wrt_key == "fb" && wrt_tag == "") {
      (*results[0]) = 2.0;
    } else {
      (*results[0]) = 0.;
    }
  }
};


class AIndependent : public EvaluatorIndependent<double> {
 public:
  using EvaluatorIndependent<double>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override {
    s.GetW<double>(my_key_, my_tag_, my_key_) = 3.0;
  }
};


class BIndependent : public EvaluatorIndependent<double> {
 public:
  using EvaluatorIndependent<double>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new BIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override {
    s.GetW<double>(my_key_, my_tag_, my_key_) = s.time(my_tag_);
  }
};


SUITE(EVALUATORS) {
  TEST(PRIMARY) {
    std::cout << "Primary Variable Test" << std::endl;
    State S;

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName("fa");

    S.Require<double>("fa", "", "fa");
    S.RequireDerivative<double>("fa", "", "fa", "");

    // Create the evaluator.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    auto fa_eval = Teuchos::rcp(new EvaluatorPrimary<double>(es_list));
    S.SetEvaluator("fa", fa_eval);

    // setup and initialize.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();
    S.GetW<double>("fa", "", "fa") = 1.0;
    S.GetRecordSetW("fa").set_initialized();
    S.Initialize();

    // provides
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", ""));     // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", "old")); // other tag
    CHECK(!S.GetEvaluator("fa").ProvidesKey("other", "")); // other key

    // dependencies -- none
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", ""));    // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", "old")); // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", "")); // not other

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // check the value and derivative
    CHECK_CLOSE(1.0, S.Get<double>("fa", ""), 1.e-10);
    CHECK_CLOSE(1.0, S.GetDerivative<double>("fa", "", "fa", ""), 1.e-10);

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fa", ""));

    // check the value and derivative are still the same
    CHECK_CLOSE(1.0, S.Get<double>("fa", ""), 1.e-10);
    CHECK_CLOSE(1.0, S.GetDerivative<double>("fa", "", "fa", ""), 1.e-10);

    // mark as changed
    auto eval = S.GetEvaluatorPtr("fa");
    auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<double>>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    // now the value is different, and so it has changed
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));

    // but the derivative is not different and has not changed
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));
  }

  TEST(SECONDARY) {
    std::cout << "Secondary Variable Test" << std::endl;
    State S;

    // make the primary.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName("fb");
    S.Require<double>("fb", "", "fb");
    auto fb_eval = Teuchos::rcp(new EvaluatorPrimary<double>(es_list));
    S.SetEvaluator("fb", fb_eval);

    // make the secondary.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    Teuchos::ParameterList ea_list;
    ea_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    ea_list.setName("fa");
    ea_list.set("tag", "");
    S.Require<double>("fa", "", "fa");
    S.RequireDerivative<double>("fa", "", "fb", "");
    auto fa_eval = Teuchos::rcp(new AEvaluator(ea_list));
    S.SetEvaluator("fa", fa_eval);

    // setup and initialize.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();

    S.GetW<double>("fa", "", "fa") = 1.0;
    S.GetW<double>("fb", "", "fb") = 3.0;

    S.GetRecordSetW("fa").set_initialized();
    S.GetRecordSetW("fb").set_initialized();
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
    CHECK_CLOSE(6.0, S.Get<double>("fa", ""), 1.e-10);
    CHECK_CLOSE(2.0, S.GetDerivative<double>("fa", "", "fb", ""), 1.e-10);

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", ""));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fb", ""));

    // check the value and derivative are still the same
    CHECK_CLOSE(6.0, S.Get<double>("fa", ""), 1.e-10);
    CHECK_CLOSE(2.0, S.GetDerivative<double>("fa", "", "fb", ""), 1.e-10);

    // change the primary and mark as changed
    S.GetW<double>("fb", "", "fb") = 14.0;
    auto eval = S.GetEvaluatorPtr("fb");
    auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<double>>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    // now the value is different, and so it has changed
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fb", ""));

    // check the values
    CHECK_CLOSE(28.0, S.Get<double>("fa", ""), 1.e-10);
    CHECK_CLOSE(2.0, S.GetDerivative<double>("fa", "", "fb", ""), 1.e-10);
  }

  TEST(INDEPENDENT_CONSTANT) {
    std::cout << "Independent Variable Test (Temporally constant)" << std::endl;
    State S;

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName("fa");
    es_list.set("constant in time", true);

    S.Require<double>("fa", "", "fa");

    // Create the evaluator.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    auto fa_eval = Teuchos::rcp(new AIndependent(es_list));
    S.SetEvaluator("fa", fa_eval);

    // setup and initialize.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();
    S.GetW<double>("fa", "", "fa") = 3.0;
    S.GetRecordSetW("fa").set_initialized();
    S.Initialize();

    // provides
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", ""));     // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", "old")); // other tag
    CHECK(!S.GetEvaluator("fa").ProvidesKey("other", "")); // other key

    // dependencies -- none
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", ""));    // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", "old")); // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", "")); // not other

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));

    // no derivatives
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // check the value and derivative
    CHECK_CLOSE(3.0, S.Get<double>("fa", ""), 1.e-10);
    CHECK_THROW(S.GetDerivative<double>("fa", "", "fa", ""), std::out_of_range);

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fa", ""));

    // check the value and derivative are still the same
    CHECK_CLOSE(3.0, S.Get<double>("fa", ""), 1.e-10);
  }

  TEST(INDEPENDENT_TEMPORALLY_CHANGING) {
    std::cout << "Independent Variable Test (Temporally varying)" << std::endl;
    State S;

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName("fa");
    es_list.set("constant in time", false);

    S.Require<double>("fa", "", "fa");

    // Create the evaluator.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    auto fa_eval = Teuchos::rcp(new BIndependent(es_list));
    S.SetEvaluator("fa", fa_eval);

    S.Require<double>("time", "", "time");

    // setup and initialize.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();

    S.GetW<double>("fa", "", "fa") = 0.0;
    S.GetRecordSetW("fa").set_initialized();
    S.Initialize();

    S.set_time(1.1);

    // provides
    CHECK(S.GetEvaluator("fa").ProvidesKey("fa", ""));     // self
    CHECK(!S.GetEvaluator("fa").ProvidesKey("fa", "old")); // other tag
    CHECK(!S.GetEvaluator("fa").ProvidesKey("other", "")); // other key

    // dependencies -- none
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", ""));    // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "fa", "old")); // not self
    CHECK(!S.GetEvaluator("fa").IsDependency(S, "other", "")); // not other

    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));

    // no derivatives
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // check the value and derivative
    CHECK_CLOSE(1.1, S.Get<double>("fa", ""), 1.e-10);
    CHECK_THROW(S.GetDerivative<double>("fa", "", "fa", ""), std::out_of_range);

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request", "fa", ""));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa").Update(S, "my_request_2"));
    CHECK(!S.GetEvaluator("fa").UpdateDerivative(S, "my_request_2", "fa", ""));

    // check the value and derivative are still the same
    CHECK_CLOSE(1.1, S.Get<double>("fa", ""), 1.e-10);

    // update time
    S.advance_time(2.0);

    // call after new time
    CHECK(S.GetEvaluator("fa").Update(S, "my_request"));
    CHECK_CLOSE(3.1, S.Get<double>("fa", ""), 1.e-10);
    CHECK(!S.GetEvaluator("fa").Update(S, "my_request"));
  }
}
