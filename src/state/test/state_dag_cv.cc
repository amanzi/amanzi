/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziComm.hh"
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
// Device Type : GPUs or CPUs??

template <class DeviceType>
class AEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  AEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.insert(std::make_pair(Key("fb"), Tag()));
    dependencies_.insert(std::make_pair(Key("fc"), Tag()));
    dependencies_.insert(std::make_pair(Key("fe"), Tag()));
    dependencies_.insert(std::make_pair(Key("fh"), Tag()));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AEvaluator(*this));
  }

  virtual std::string getType() const override { return "AEvaluator"; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fb_c = S.Get<CompositeVector>("fb").viewComponent("cell", false);
    const auto fc_c = S.Get<CompositeVector>("fc").viewComponent("cell", false);
    const auto fe_c = S.Get<CompositeVector>("fe").viewComponent("cell", false);
    const auto fh_c = S.Get<CompositeVector>("fh").viewComponent("cell", false);
    Kokkos::parallel_for(
      1, KOKKOS_LAMBDA(const int i) {
        result_c(i, 0) = 2 * fb_c(i, 0) + fc_c(i, 0) * fe_c(i, 0) * fh_c(i, 0);
      });
  }

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fc_c = S.Get<CompositeVector>("fc").viewComponent("cell");
    const auto fe_c = S.Get<CompositeVector>("fe").viewComponent("cell");
    const auto fh_c = S.Get<CompositeVector>("fh").viewComponent("cell");

    if (wrt_key == "fb") {
      Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2.0; });
    } else if (wrt_key == "fc") {
      Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) { result_c(i, 0) = fe_c(i, 0) * fh_c(i, 0); });
    } else if (wrt_key == "fe") {
      Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) { result_c(i, 0) = fc_c(i, 0) * fh_c(i, 0); });
    } else if (wrt_key == "fh") {
      Kokkos::parallel_for(
        1, KOKKOS_LAMBDA(const int i) { result_c(i, 0) = fc_c(i, 0) * fe_c(i, 0); });
    }
  }
};

/* ******************************************************************
 * Equation C = 2*D + G
 ****************************************************************** */
class CEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  CEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.insert(std::make_pair(Key("fd"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("fg"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new CEvaluator(*this));
  }
  virtual std::string getType() const override { return "CEvaluator"; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fd_c = S.Get<CompositeVector>("fd").viewComponent("cell", false);
    const auto fg_c = S.Get<CompositeVector>("fg").viewComponent("cell", false);
    Kokkos::parallel_for(
      result_c.extent(0),
      KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2 * fd_c(i, 0) + fg_c(i, 0); });
  }

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    if (wrt_key == "fd") {
      Kokkos::parallel_for(result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2.; });
    } else if (wrt_key == "fg") {
      Kokkos::parallel_for(result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 1.; });
    }
  }
};

/* ******************************************************************
 * Equation D = 2*G
 ****************************************************************** */
class DEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  DEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.insert(std::make_pair(Key("fg"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new DEvaluator(*this));
  }
  virtual std::string getType() const override { return "DEvaluator"; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fg_c = S.Get<CompositeVector>("fg").viewComponent("cell", false);
    Kokkos::parallel_for(
      result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2 * fg_c(i, 0); });
  }

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    if (wrt_key == "fg") {
      Kokkos::parallel_for(result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2.; });
    }
  }
};

/* ******************************************************************
 * Equation E = D*F
 ****************************************************************** */
class EEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  EEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.insert(std::make_pair(Key("fd"), Tags::DEFAULT));
    dependencies_.insert(std::make_pair(Key("ff"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EEvaluator(*this));
  }
  virtual std::string getType() const override { return "EEvaluator"; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fd_c = S.Get<CompositeVector>("fd").viewComponent("cell", false);
    const auto ff_c = S.Get<CompositeVector>("ff").viewComponent("cell", false);
    Kokkos::parallel_for(
      result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = fd_c(i, 0) * ff_c(i, 0); });
  }

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fd_c = S.Get<CompositeVector>("fd").viewComponent("cell", false);
    const auto ff_c = S.Get<CompositeVector>("ff").viewComponent("cell", false);

    if (wrt_key == "fd") {
      Kokkos::parallel_for(
        result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = ff_c(i, 0); });
    } else if (wrt_key == "ff") {
      Kokkos::parallel_for(
        result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = fd_c(i, 0); });
    }
  }
};

/* ******************************************************************
 * Equation F = 2*G
 ****************************************************************** */
class FEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  FEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.insert(std::make_pair(Key("fg"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new FEvaluator(*this));
  }
  virtual std::string getType() const override { return "FEvaluator"; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto fg_c = S.Get<CompositeVector>("fg").viewComponent("cell", false);
    Kokkos::parallel_for(
      result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2 * fg_c(i, 0); });
  }

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    if (wrt_key == "fg") {
      Kokkos::parallel_for(result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2.; });
    }
  }
};

/* ******************************************************************
 * Equation H = 2*F
 ****************************************************************** */
class HEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  HEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
  {
    dependencies_.insert(std::make_pair(Key("ff"), Tags::DEFAULT));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new HEvaluator(*this));
  }
  virtual std::string getType() const override { return "HEvaluator"; }

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    const auto ff_c = S.Get<CompositeVector>("ff").viewComponent("cell", false);
    Kokkos::parallel_for(
      result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2. * ff_c(i, 0); });
  }

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    auto result_c = results[0]->viewComponent("cell", false);
    if (wrt_key == "ff") {
      Kokkos::parallel_for(result_c.extent(0), KOKKOS_LAMBDA(const int i) { result_c(i, 0) = 2.; });
    }
  }
};

class make_state {
 public:
  make_state()
  {
    auto es_list = Teuchos::rcp(new Teuchos::ParameterList());
    es_list->sublist("verbose object").set<std::string>("verbosity level", "extreme");
    es_list->set("tag", "");
    auto ep_list = Teuchos::rcp(new Teuchos::ParameterList());
    ep_list->sublist("verbose object").set<std::string>("verbosity level", "extreme");

    // create a mesh
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // create a state
    // State S;
    S.RegisterDomainMesh(mesh);

    // Primary fields
    ep_list->setName("fb");
    // -- field B and its evaluator
    S.Require<CompositeVector, CompositeVectorSpace>("fb", Tags::DEFAULT, "fb")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fb", Tags::DEFAULT, fb_eval);

    // -- field G and its evaluator
    ep_list->setName("fg");
    S.Require<CompositeVector, CompositeVectorSpace>("fg", Tags::DEFAULT, "fg")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    fg_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fg", Tags::DEFAULT, fg_eval);

    // Secondary fields
    // --  D and its evaluator
    es_list->setName("fd");
    S.Require<CompositeVector, CompositeVectorSpace>("fd", Tags::DEFAULT, "fd")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    fd_eval = Teuchos::rcp(new DEvaluator(es_list));
    S.SetEvaluator("fd", Tags::DEFAULT, fd_eval);

    // --  C and its evaluator
    es_list->setName("fc");
    S.Require<CompositeVector, CompositeVectorSpace>("fc", Tags::DEFAULT, "fc")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    fc_eval = Teuchos::rcp(new CEvaluator(es_list));
    S.SetEvaluator("fc", Tags::DEFAULT, fc_eval);

    // --  F and its evaluator
    es_list->setName("ff");
    S.Require<CompositeVector, CompositeVectorSpace>("ff", Tags::DEFAULT, "ff")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    ff_eval = Teuchos::rcp(new FEvaluator(es_list));
    S.SetEvaluator("ff", Tags::DEFAULT, ff_eval);

    // --  E and its evaluator
    es_list->setName("fe");
    S.Require<CompositeVector, CompositeVectorSpace>("fe", Tags::DEFAULT, "fe")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fe", Tags::DEFAULT, "fg", Tags::DEFAULT);
    fe_eval = Teuchos::rcp(new EEvaluator(es_list));
    S.SetEvaluator("fe", Tags::DEFAULT, fe_eval);

    // --  H and its evaluator
    es_list->setName("fh");
    S.Require<CompositeVector, CompositeVectorSpace>("fh", Tags::DEFAULT, "fh")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    fh_eval = Teuchos::rcp(new HEvaluator(es_list));
    S.SetEvaluator("fh", Tags::DEFAULT, fh_eval);

    // --  A and its evaluator
    es_list->setName("fa");
    S.Require<CompositeVector, CompositeVectorSpace>("fa", Tags::DEFAULT, "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", Tags::DEFAULT, "fb", Tags::DEFAULT);
    S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
      "fa", Tags::DEFAULT, "fg", Tags::DEFAULT);
    fa_eval = Teuchos::rcp(new AEvaluator<DefaultDevice>(es_list));
    S.SetEvaluator("fa", Tags::DEFAULT, fa_eval);

    // Setup fields initialize
    S.Setup();
    S.GetW<CompositeVector>("fb", "fb").putScalar(2.0);
    S.GetRecordW("fb", "fb").set_initialized();
    S.GetW<CompositeVector>("fg", "fg").putScalar(3.0);
    S.GetRecordW("fg", "fg").set_initialized();
    S.Initialize();
  }

 public:
  State S;
  Teuchos::RCP<AEvaluator<DefaultDevice>> fa_eval;
  Teuchos::RCP<CEvaluator> fc_eval;
  Teuchos::RCP<DEvaluator> fd_eval;
  Teuchos::RCP<EEvaluator> fe_eval;
  Teuchos::RCP<FEvaluator> ff_eval;
  Teuchos::RCP<HEvaluator> fh_eval;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> fb_eval, fg_eval;
};

SUITE(DAG)
{
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS)
  {
    // check initialized properly
    CHECK_CLOSE(
      2.0, (S.Get<CompositeVector>("fb").viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0), 1e-12);
    CHECK_CLOSE(
      3.0, (S.Get<CompositeVector>("fg").viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0), 1e-12);

    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    CHECK_CLOSE(
      6484.0, (S.Get<CompositeVector>("fa").viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0), 1e-12);
    CHECK(changed);


    // check intermediate steps got updated too
    CHECK_CLOSE(
      6.0, (S.Get<CompositeVector>("fd").viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0), 1e-12);

    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fb", Tags::DEFAULT);
    CHECK_CLOSE(2.0,
                (S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fb", Tags::DEFAULT)
                   .viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0),
                1e-12);
    CHECK(changed);

    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", Tags::DEFAULT);
    CHECK_CLOSE(8640.0,
                (S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT)
                   .viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0),
                1e-12);
    CHECK(changed);

    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    changed = fe_eval->UpdateDerivative(S, "fe", "fg", Tags::DEFAULT);
    CHECK_CLOSE(24.0,
                (S.GetDerivative<CompositeVector>("fe", Tags::DEFAULT, "fg", Tags::DEFAULT)
                   .viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0),
                1e-12);
    CHECK(changed);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", Tags::DEFAULT);
    CHECK_CLOSE(8640.0,
                (S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT)
                   .viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0),
                1e-12);
    CHECK(!changed);

    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetChanged();
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", Tags::DEFAULT);
    CHECK_CLOSE(8640.0,
                (S.GetDerivative<CompositeVector>("fa", Tags::DEFAULT, "fg", Tags::DEFAULT)
                   .viewComponent<MemSpace_kind::HOST>("cell", false))(0, 0),
                1e-12);
    CHECK(changed);
  }
}
