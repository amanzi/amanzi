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
class AEvaluator : public EvaluatorSecondary<double,NullFactory> {
 public:
  AEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<double,NullFactory>(plist) {
    my_key_ = "fa";
    dependencies_.insert("fb");
    dependencies_.insert("fc");
    dependencies_.insert("fe");
    dependencies_.insert("fh");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this)); }

  virtual void Evaluate_(State& S,
          double& result) override {

    auto& fb = S.Get<double>("fb");
    auto& fc = S.Get<double>("fc");
    auto& fe = S.Get<double>("fe");
    auto& fh = S.Get<double>("fh");
    result = 2 * fb + fc * fe * fh;
  }

  virtual void EvaluatePartialDerivative_(State& S,
          const Key& wrt_key, double& result) override {
    auto& fc = S.Get<double>("fc");
    auto& fe = S.Get<double>("fe");
    auto& fh = S.Get<double>("fh");
    
    if (wrt_key == "fb") {
      result = 2.0;
    } else if (wrt_key == "fc") {
      result = fe * fh;
    } else if (wrt_key == "fe") {
      result = fc * fh;
    } else if (wrt_key == "fh") {
      result = fc * fe;
    }
  }
};


/* ******************************************************************
* Equation C = 2*D + G
****************************************************************** */
class CEvaluator : public EvaluatorSecondary<double,NullFactory> {
 public:
  CEvaluator(Teuchos::ParameterList& plist) : 
    EvaluatorSecondary<double,NullFactory>(plist) {
    my_key_ = "fc";
    dependencies_.insert("fd");
    dependencies_.insert("fg");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new CEvaluator(*this)); }

  virtual void Evaluate_(State& S,
          double& result) override {
    auto& fd = S.Get<double>("fd");
    auto& fg = S.Get<double>("fg");
    result = 2 * fd + fg;
  }

  virtual void EvaluatePartialDerivative_(State& S,
          const Key& wrt_key, double& result) override {

    if (wrt_key == "fd") {
      result = 2.;
    } else if (wrt_key == "fg") {
      result = 1.;
    }
  }
};


/* ******************************************************************
* Equation D = 2*G
****************************************************************** */
class DEvaluator : public EvaluatorSecondary<double,NullFactory> {
 public:
  DEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<double,NullFactory>(plist) {
    my_key_ = "fd";
    dependencies_.insert("fg");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DEvaluator(*this)); }

  virtual void Evaluate_(State& S,
          double& result) override {
    auto& fg = S.Get<double>("fg");
    result = 2*fg;
  }

  virtual void EvaluatePartialDerivative_(State& S,
          const Key& wrt_key, double& result) override {
    if (wrt_key == "fg") {
      result = 2.;
    }
  }
};


/* ******************************************************************
* Equation E = D*F
****************************************************************** */
class EEvaluator : public EvaluatorSecondary<double,NullFactory> {
 public:
  EEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<double,NullFactory>(plist) {
    my_key_ = "fe";
    dependencies_.insert("fd");
    dependencies_.insert("ff");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EEvaluator(*this)); }

  virtual void Evaluate_(State& S,
          double& result) override {
    auto& fd = S.Get<double>("fd");
    auto& ff = S.Get<double>("ff");
    result = fd * ff;
  }

  virtual void EvaluatePartialDerivative_(State& S,
          const Key& wrt_key, double& result) override {
    auto& fd = S.Get<double>("fd");
    auto& ff = S.Get<double>("ff");

    if (wrt_key == "fd") {
      result = ff;
    } else if (wrt_key == "ff") {
      result = fd;
    }
  }
};


/* ******************************************************************
* Equation F = 2*G
****************************************************************** */
class FEvaluator : public EvaluatorSecondary<double,NullFactory> {
 public:
  FEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<double,NullFactory>(plist) {
    my_key_ = "ff";
    dependencies_.insert("fg");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new FEvaluator(*this)); }

  virtual void Evaluate_(State& S,
          double& result) override {
    auto& fg = S.Get<double>("fg");
    result = 2*fg;
  }

  virtual void EvaluatePartialDerivative_(State& S,
          const Key& wrt_key, double& result) override {

    if (wrt_key == "fg") {
      result = 2.;
    }
  }
};


/* ******************************************************************
* Equation H = 2*F
****************************************************************** */
class HEvaluator : public EvaluatorSecondary<double,NullFactory> {
 public:
  HEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondary<double,NullFactory>(plist) {
    my_key_ = "fh";
    dependencies_.insert("ff");
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new HEvaluator(*this)); }

  virtual void Evaluate_(State& S,
          double& result) override {
    auto& ff = S.Get<double>("ff");
    result = 2.*ff;
  }

  virtual void EvaluatePartialDerivative_(State& S,
          const Key& wrt_key, double& result) override {
    if (wrt_key == "ff") {
      result = 2.;
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
    S.Require<double>("fa", "", "fa");
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S.SetEvaluator("fa", fa_eval);

    // --  C and its evaluator
    es_list.setName("fc");
    S.Require<double>("fc", "", "fc");
    fc_eval = Teuchos::rcp(new CEvaluator(es_list));
    S.SetEvaluator("fc", fc_eval);

    // --  D and its evaluator
    es_list.setName("fd");
    S.Require<double>("fd", "", "fd");
    fd_eval = Teuchos::rcp(new DEvaluator(es_list));
    S.SetEvaluator("fd", fd_eval);

    // --  E and its evaluator
    es_list.setName("fe");
    S.Require<double>("fe", "", "fe");
    fe_eval = Teuchos::rcp(new EEvaluator(es_list));
    S.SetEvaluator("fe", fe_eval);

    // --  F and its evaluator
    es_list.setName("ff");
    S.Require<double>("ff", "", "ff");
    ff_eval = Teuchos::rcp(new FEvaluator(es_list));
    S.SetEvaluator("ff", ff_eval);

    // --  H and its evaluator
    es_list.setName("fh");
    S.Require<double>("fh", "", "fh");
    fh_eval = Teuchos::rcp(new HEvaluator(es_list));
    S.SetEvaluator("fh", fh_eval);

    // Primary fields
    ep_list.setName("fb");
    // -- field B and its evaluator
    S.Require<double>("fb", "", "fb");
    fb_eval = Teuchos::rcp(new EvaluatorPrimary(ep_list));
    S.SetEvaluator("fb", fb_eval);

    // -- field G and its evaluator
    ep_list.setName("fg");
    S.Require<double>("fg", "", "fg");
    fg_eval = Teuchos::rcp(new EvaluatorPrimary(ep_list));
    S.SetEvaluator("fg", fg_eval);

    // Setup fields initialize
    S.Setup();
    S.GetW<double>("fb","fb") = 2.0;
    S.GetRecordW("fb", "fb").set_initialized();
    S.GetW<double>("fg","fg") = 3.0;
    S.GetRecordW("fg", "fg").set_initialized();
    S.Initialize();
  }

 public:
  State S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<CEvaluator> fc_eval;
  Teuchos::RCP<DEvaluator> fd_eval;
  Teuchos::RCP<EEvaluator> fe_eval;
  Teuchos::RCP<FEvaluator> ff_eval;
  Teuchos::RCP<HEvaluator> fh_eval;
  Teuchos::RCP<EvaluatorPrimary> fb_eval, fg_eval;
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
    changed = fa_eval->UpdateDerivative(S, "fa", "fb");
    CHECK_CLOSE(2.0, S.Get<double>("dfa_dfb"), 1e-12);
    CHECK(changed);

    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg");
    CHECK_CLOSE(8640.0, S.Get<double>("dfa_dfg"), 1e-12);
    CHECK(changed);

    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    changed = fe_eval->UpdateDerivative(S, "fe", "fg");
    CHECK_CLOSE(24.0, S.Get<double>("dfe_dfg"), 1e-12);
    CHECK(changed);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg");
    CHECK_CLOSE(8640.0, S.Get<double>("dfa_dfg"), 1e-12);
    CHECK(!changed);

    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetChanged();
    changed = fa_eval->UpdateDerivative(S, "fa", "fg");
    CHECK_CLOSE(8640.0, S.Get<double>("dfa_dfg"), 1e-12);
    CHECK(changed);
  }
}

