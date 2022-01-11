/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "IO.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "EvaluatorPrimary.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

/*
  We will build the following dependencies tree:
    A -> {B, C, E, H}
    C -> {D, G}
    E -> {D, F}
    H -> F
    D -> G
    F -> G

  Primary fields are B=2 and G=3. The equations are
    A = 2*B + C*E*H = 6484
    C = 2*D + G     = 15
    E = D*F         = 36
    H = 2*F         = 12
    D = 2*G         = 6
    F = 2*G         = 6

  Derivatives are
    dA/dB = 2
    dA/dG = 8640

  WARNING: derivative of secondary field wrt to secondary field is
  not well defined. The code may throw an exception since
  intermediate derivatives are not saved.
*/

/* ******************************************************************
 * Equation A = 2*B + C*E*H
 ****************************************************************** */
class AEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  AEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.insert(std::make_pair(Key("fb"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("fc"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("fe"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("fh"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<double*>& results) override {

    auto& fb = S.Get<double>("fb");
    auto& fc = S.Get<double>("fc");
    auto& fe = S.Get<double>("fe");
    auto& fh = S.Get<double>("fh");
    (*results[0]) = 2 * fb + fc * fe * fh;
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<double*>& results) override {
    auto& fc = S.Get<double>("fc");
    auto& fe = S.Get<double>("fe");
    auto& fh = S.Get<double>("fh");

    if (wrt_key == "fb") {
      (*results[0]) = 2.0;
    } else if (wrt_key == "fc") {
      (*results[0]) = fe * fh;
    } else if (wrt_key == "fe") {
      (*results[0]) = fc * fh;
    } else if (wrt_key == "fh") {
      (*results[0]) = fc * fe;
    }
  }
};

/* ******************************************************************
 * Equation C = 2*D + G
 ****************************************************************** */
class CEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  CEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.insert(std::make_pair(Key("fd"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("fg"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new CEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<double*>& results) override {
    auto& fd = S.Get<double>("fd");
    auto& fg = S.Get<double>("fg");
    (*results[0]) = 2 * fd + fg;
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<double*>& results) override {

    if (wrt_key == "fd") {
      (*results[0]) = 2.;
    } else if (wrt_key == "fg") {
      (*results[0]) = 1.;
    }
  }
};

/* ******************************************************************
 * Equation D = 2*G
 ****************************************************************** */
class DEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  DEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.insert(std::make_pair(Key("fg"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<double*>& results) override {
    auto& fg = S.Get<double>("fg");
    (*results[0]) = 2 * fg;
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<double*>& results) override {
    if (wrt_key == "fg") {
      (*results[0]) = 2.;
    }
  }
};

/* ******************************************************************
 * Equation E = D*F
 ****************************************************************** */
class EEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  EEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.insert(std::make_pair(Key("fd"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("ff"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<double*>& results) override {
    auto& fd = S.Get<double>("fd");
    auto& ff = S.Get<double>("ff");
    (*results[0]) = fd * ff;
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<double*>& results) override {
    auto& fd = S.Get<double>("fd");
    auto& ff = S.Get<double>("ff");

    if (wrt_key == "fd") {
      (*results[0]) = ff;
    } else if (wrt_key == "ff") {
      (*results[0]) = fd;
    }
  }
};

/* ******************************************************************
 * Equation F = 2*G
 ****************************************************************** */
class FEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  FEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.insert(std::make_pair(Key("fg"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new FEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<double*>& results) override {
    auto& fg = S.Get<double>("fg");
    (*results[0]) = 2 * fg;
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<double*>& results) override {

    if (wrt_key == "fg") {
      (*results[0]) = 2.;
    }
  }
};

/* ******************************************************************
 * Equation H = 2*F
 ****************************************************************** */
class HEvaluator : public EvaluatorSecondaryMonotype<double> {
public:
  HEvaluator(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<double>(plist) {
    dependencies_.insert(std::make_pair(Key("ff"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new HEvaluator(*this));
  }

  virtual void Evaluate_(const State& S, const std::vector<double*>& results) override {
    auto& ff = S.Get<double>("ff");
    (*results[0]) = 2. * ff;
  }

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<double*>& results) override {
    if (wrt_key == "ff") {
      (*results[0]) = 2.;
    }
  }
};

class make_state {
public:
  make_state() {
    Teuchos::ParameterList es_list, ep_list;
    es_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    ep_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");

    // Secondary fields
    // --  A and its evaluator
    es_list.setName("fa");
    es_list.set("tag", "");
    S.Require<double>("fa", Tags::DEFAULT, "fa");
    S.RequireDerivative<double>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    S.RequireDerivative<double>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT);
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S.SetEvaluator("fa", fa_eval);

    // --  C and its evaluator
    es_list.setName("fc");
    S.Require<double>("fc", Tags::DEFAULT, "fc");
    fc_eval = Teuchos::rcp(new CEvaluator(es_list));
    S.SetEvaluator("fc", fc_eval);

    // --  D and its evaluator
    es_list.setName("fd");
    S.Require<double>("fd", Tags::DEFAULT, "fd");
    fd_eval = Teuchos::rcp(new DEvaluator(es_list));
    S.SetEvaluator("fd", fd_eval);

    // --  E and its evaluator
    es_list.setName("fe");
    S.Require<double>("fe", Tags::DEFAULT, "fe");
    S.RequireDerivative<double>("fe", Tags::DEFAULT, "fg", Tags::DEFAULT);
    fe_eval = Teuchos::rcp(new EEvaluator(es_list));
    S.SetEvaluator("fe", fe_eval);

    // --  F and its evaluator
    es_list.setName("ff");
    S.Require<double>("ff", Tags::DEFAULT, "ff");
    ff_eval = Teuchos::rcp(new FEvaluator(es_list));
    S.SetEvaluator("ff", ff_eval);

    // --  H and its evaluator
    es_list.setName("fh");
    S.Require<double>("fh", Tags::DEFAULT, "fh");
    fh_eval = Teuchos::rcp(new HEvaluator(es_list));
    S.SetEvaluator("fh", fh_eval);

    // Primary fields
    ep_list.setName("fb");
    // -- field B and its evaluator
    S.Require<double>("fb", Tags::DEFAULT, "fb");
    fb_eval = Teuchos::rcp(new EvaluatorPrimary<double>(ep_list));
    S.SetEvaluator("fb", fb_eval);

    // -- field G and its evaluator
    ep_list.setName("fg");
    S.Require<double>("fg", Tags::DEFAULT, "fg");
    fg_eval = Teuchos::rcp(new EvaluatorPrimary<double>(ep_list));
    S.SetEvaluator("fg", fg_eval);

    // Setup fields initialize
    S.Setup();
    S.GetW<double>("fb", "fb") = 2.0;
    S.GetRecordSetW("fb").initializeTags();
    S.GetW<double>("fg", "fg") = 3.0;
    S.GetRecordSetW("fg").initializeTags();
    S.Initialize();

    Teuchos::ParameterList plist;
    plist.sublist("verbose object").set("verbosity level", "extreme");

    VerboseObject vo("State", plist);
    WriteStateStatistics(S, vo);
  }

public:
  State S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<CEvaluator> fc_eval;
  Teuchos::RCP<DEvaluator> fd_eval;
  Teuchos::RCP<EEvaluator> fe_eval;
  Teuchos::RCP<FEvaluator> ff_eval;
  Teuchos::RCP<HEvaluator> fh_eval;
  Teuchos::RCP<EvaluatorPrimary<double>> fb_eval, fg_eval;
};


SUITE(DAG) {
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS) {
    // check initialized properly
    CHECK_CLOSE(2.0, S.Get<double>("fb"), 1e-12);
    CHECK_CLOSE(3.0, S.Get<double>("fg"), 1e-12);

    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    CHECK_CLOSE(6484.0, S.Get<double>("fa"), 1e-12);
    CHECK(changed);

    // check intermediate steps got updated too
    CHECK_CLOSE(6.0, S.Get<double>("fd"), 1e-12);

    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fb", Tags::DEFAULT);
    CHECK_CLOSE(2.0, S.GetDerivative<double>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT), 1e-12);
    CHECK(changed);

    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", Tags::DEFAULT);
    CHECK_CLOSE(8640.0, S.GetDerivative<double>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT), 1e-12);
    CHECK(changed);

    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    changed = fe_eval->UpdateDerivative(S, "fe", "fg", Tags::DEFAULT);
    CHECK_CLOSE(24.0, S.GetDerivative<double>("fe", Tags::DEFAULT, "fg", Tags::DEFAULT), 1e-12);
    CHECK(changed);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", Tags::DEFAULT);
    CHECK_CLOSE(8640.0, S.GetDerivative<double>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT), 1e-12);
    CHECK(!changed);

    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetChanged();
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", Tags::DEFAULT);
    CHECK_CLOSE(8640.0, S.GetDerivative<double>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT), 1e-12);
    CHECK(changed);
  }
}
