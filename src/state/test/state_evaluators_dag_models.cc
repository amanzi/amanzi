/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/


//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "dag_models.hh"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "evaluator/EvaluatorModel.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;


class make_state {
 public:
  make_state()
  {
    // create a mesh
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // create a state
    // State S;
    S.RegisterDomainMesh(mesh);

    // Secondary fields
    fa_eval = requireSecondary<AModel>("A", { "B", "G" });
    fc_eval = requireSecondary<CModel>("C", {});
    fd_eval = requireSecondary<DModel>("D", {});
    fe_eval = requireSecondary<EModel>("E", { "G", });
    ff_eval = requireSecondary<FModel>("F", {});
    fh_eval = requireSecondary<HModel>("H", {});

    // Primary fields
    fb_eval = requirePrimary("B");
    fg_eval = requirePrimary("G");

    // Setup fields initialize
    S.Setup();
    S.GetW<CompositeVector>("B", "B").putScalar(2.0);
    S.GetRecordW("B", "B").set_initialized();
    S.GetW<CompositeVector>("G", "G").putScalar(3.0);
    S.GetRecordW("G", "G").set_initialized();
    S.Initialize();
  }

  template<template<class,class> class Model>
  Teuchos::RCP<Evaluator>
  requireSecondary(const std::string& name, const std::vector<std::string>& derivs) {
    Teuchos::ParameterList es_list(name);
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName(name);
    es_list.set("tag", "");
    S.Require<CompositeVector,CompositeVectorSpace>(name, "", name)
        .SetMesh(S.GetMesh())
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    for (const auto& deriv : derivs)
      S.RequireDerivative<CompositeVector,CompositeVectorSpace>(name, "", deriv, "");
    
    //  S.RequireDerivative<double>("fa", "", "fb", "");
    //  S.RequireDerivative<double>("fa", "", "fg", "");
    auto f_eval = Teuchos::rcp(new EvaluatorModel_CompositeVector<Model>(es_list));
    S.SetEvaluator(name, f_eval);
    return f_eval;
  }


  Teuchos::RCP<EvaluatorPrimary_>
  requirePrimary(const std::string& name) {
    Teuchos::ParameterList es_list(name);
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName(name);
    es_list.set("tag", "");
    S.Require<CompositeVector,CompositeVectorSpace>(name, "", name)
        .SetMesh(S.GetMesh())
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    
    //  S.RequireDerivative<double>("fa", "", "fb", "");
    //  S.RequireDerivative<double>("fa", "", "fg", "");
    auto f_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(es_list));
    S.SetEvaluator(name, f_eval);
    return f_eval;
  }

  void check_close(double val, const std::string& name) {
    auto cvv = S.Get<CompositeVector>(name, "").ViewComponent<AmanziDefaultHost>("cell", 0, false);
    CHECK_CLOSE(val, cvv(0), 1.e-10);
  }

  void check_close_deriv(double val, const std::string& name, const std::string& wrt) {
    auto cvv = S.GetDerivative<CompositeVector>(name, "", wrt, "").ViewComponent<AmanziDefaultHost>("cell", 0, false);
    CHECK_CLOSE(val, cvv(0), 1.e-10);
  }
  
  
 public:
  State S;
  Teuchos::RCP<Evaluator> fa_eval;
  Teuchos::RCP<Evaluator> fc_eval;
  Teuchos::RCP<Evaluator> fd_eval;
  Teuchos::RCP<Evaluator> fe_eval;
  Teuchos::RCP<Evaluator> ff_eval;
  Teuchos::RCP<Evaluator> fh_eval;
  Teuchos::RCP<EvaluatorPrimary_> fb_eval, fg_eval;
};

SUITE(DAG)
{
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS)
  {
    // check initialized properly
    check_close(2.0,"B");
    check_close(3.0,"G");

    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    check_close(6484.0,"A");
    CHECK(changed);

    // check intermediate steps got updated too
    check_close(6.0, "D");

    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "A", "B", "");
    check_close_deriv(2.0, "A", "B");
    CHECK(changed);

    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "A", "G", "");
    check_close_deriv(8640.0, "A", "G");
    CHECK(changed);

    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    changed = fe_eval->UpdateDerivative(S, "E", "G", "");
    check_close_deriv(24.0, "E", "G");
    CHECK(changed);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "A", "G", "");
    check_close_deriv(8640.0, "A", "G");
    CHECK(!changed);

    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetChanged();
    changed = fa_eval->UpdateDerivative(S, "A", "G", "");
    check_close_deriv(8640.0, "A", "G");
    CHECK(changed);
  }
}
