/*
  State

  Authors: Ethan Coon
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "evaluator/EvaluatorSecondary.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/* ******************************************************************
* Equation A = 2*B
****************************************************************** */
class AEvaluator : public EvaluatorSecondary<double> {
 public:
  AEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<double>(plist) {
    my_key_ = std::string("fa");
    dependencies_.insert(std::string("fb"));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this)); };

  virtual void Evaluate_(const State& S,
          double& result) override {
    auto& fb = S.Get<double>("fb");
    result = 2 * fb;
  }

  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, double& result) override {
    if (wrt_key == "fb") {
      result = 2.0;
    }
  }
};


SUITE(EVALUATORS) {
  TEST(PRIMARY) {
    State S;

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    es_list.setName("fa");

    S.Require<double>("fa", "", "fa");

    auto fa_eval = Teuchos::rcp(new EvaluatorPrimary(es_list));
    S.SetEvaluator("fa", fa_eval);

    // Setup fields and marked as initialized.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();
    S.GetW<double>("fa", "", "fa") = 1.0;
    S.GetRecordW("fa", "fa").set_initialized();
    S.Initialize();

    // provides
    CHECK(S.GetEvaluator("fa")->ProvidesKey("fa")); // self
    CHECK(!S.GetEvaluator("fa")->ProvidesKey("other")); // self
    
    // dependencies -- none
    CHECK(!S.GetEvaluator("fa")->IsDependency(S, "fa")); // not self
    CHECK(!S.GetEvaluator("fa")->IsDependency(S, "other")); // not other
    
    // check first call is always "changed"
    CHECK(S.GetEvaluator("fa")->Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa")->UpdateDerivative(S, "my_request", "fa"));

    // second call should not be changed
    CHECK(!S.GetEvaluator("fa")->Update(S, "my_request"));
    CHECK(!S.GetEvaluator("fa")->UpdateDerivative(S, "my_request", "fa"));

    // but first call with new request should be
    CHECK(S.GetEvaluator("fa")->Update(S, "my_request_2"));
    CHECK(S.GetEvaluator("fa")->UpdateDerivative(S, "my_request_2", "fa"));

    // mark as changed
    auto eval = S.GetEvaluator("fa");
    auto eval_p = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(eval);
    CHECK(eval_p.get());
    eval_p->SetChanged();

    CHECK(S.GetEvaluator("fa")->Update(S, "my_request"));
    CHECK(S.GetEvaluator("fa")->UpdateDerivative(S, "my_request", "fa"));
  }


}

