/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "Executor.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;



template<class DeviceType>
class AModel {
 public:
  AModel(OutputVector_type<DeviceType> A, InputVector_type<DeviceType> B, InputVector_type<DeviceType> C, 
	 InputVector_type<DeviceType> E, InputVector_type<DeviceType> H, Teuchos::ParameterList &plist) :
    A_(A), B_(B), C_(C), E_(E), H_(H) {}

  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    A_(i) = 2 * B_(i)  + C_(i) * E_(i) * H_(i);
  }

  class dAdB {};
  KOKKOS_INLINE_FUNCTION void operator() (dAdB, const int i) const {
    A_(i) = 2.0;
  }

  class dAdC{};
  KOKKOS_INLINE_FUNCTION void operator() (dAdC, const int i) const {
    A_(i) = E_(i) * H_(i);
  }
  class dAdE{};
  KOKKOS_INLINE_FUNCTION void operator() (dAdE, const int i) const {
    A_(i) = C_(i) * H_(i);
  } 
  
  class dAdH{};
  KOKKOS_INLINE_FUNCTION void operator() (dAdH, const int i) const {
    A_(i) = C_(i) * E_(i);
  }

private:
  OutputVector_type<DeviceType> A_;
  InputVector_type<DeviceType> B_;
  InputVector_type<DeviceType> C_;
  InputVector_type<DeviceType> E_;
  InputVector_type<DeviceType> H_;
};



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


class AEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
public:
  AEvaluator(Teuchos::ParameterList &plist)
    : EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace>(plist) {
    dependencies_.emplace_back(std::make_pair(Key("fb"), Key()));
    dependencies_.emplace_back(std::make_pair(Key("fc"), Key()));
    dependencies_.emplace_back(std::make_pair(Key("fe"), Key()));
    dependencies_.emplace_back(std::make_pair(Key("fh"), Key()));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this));
  }

  virtual void Evaluate_(const State &S, const std::vector<CompositeVector*> &results) override {
    auto result_c = results[0]->ViewComponent("cell",0,false);
    auto fb_c = S.Get<CompositeVector>("fb").ViewComponent("cell",0,false);
    auto fc_c = S.Get<CompositeVector>("fc").ViewComponent("cell",0,false);
    auto fe_c = S.Get<CompositeVector>("fe").ViewComponent("cell",0,false);
    auto fh_c = S.Get<CompositeVector>("fh").ViewComponent("cell",0,false);

    AModel<AmanziDefaultDevice> model(result_c, fb_c, fc_c, fe_c, fh_c, plist_);
    ExecuteModel("A", model, 0, result_c.extent(0));
  }

  virtual void EvaluatePartialDerivative_(const State &S, 
					  const Key &wrt_key,
                                          const Key &wrt_tag,
                                          const std::vector<CompositeVector*> &results) override {
    auto result = results[0]->ViewComponent("cell",0,false);
    auto fb = S.Get<CompositeVector>("fb").ViewComponent("cell",0,false);
    auto fc = S.Get<CompositeVector>("fc").ViewComponent("cell",0,false);
    auto fe = S.Get<CompositeVector>("fe").ViewComponent("cell",0,false);
    auto fh = S.Get<CompositeVector>("fh").ViewComponent("cell",0,false);
    AModel<AmanziDefaultDevice> model(result, fb, fc, fe, fh, plist_);
    if (wrt_key == "fb") {
      ExecuteModel("dAdB",model, AModel<AmanziDefaultDevice>::dAdB(), 0,result.extent(0));
    } else if (wrt_key == "fc") {
      ExecuteModel("dAdC",model, AModel<AmanziDefaultDevice>::dAdC(), 0,result.extent(0));
    } else if (wrt_key == "fe") {
      ExecuteModel("dAdE",model, AModel<AmanziDefaultDevice>::dAdE(), 0,result.extent(0));
    } else if (wrt_key == "fh") {
      ExecuteModel("dAdH",model, AModel<AmanziDefaultDevice>::dAdH(), 0,result.extent(0));
    }
  }
};


class make_state {
public:
  make_state() {
    Teuchos::ParameterList es_list, ep_list;
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    ep_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");

    // create a mesh
    auto comm = Comm_ptr_type( new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    MeshFactory meshfac(comm);
    auto mesh = meshfac(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // create a state
    S.RegisterDomainMesh(mesh);

    // Secondary fields
    // --  A and its evaluator
    es_list.setName("fa");
    es_list.set("tag", "");
    S.Require<CompositeVector,CompositeVectorSpace>("fa", "", "fa")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S.RequireDerivative<CompositeVector,CompositeVectorSpace>("fa", "", "fb", "");
    S.RequireDerivative<CompositeVector,CompositeVectorSpace>("fa", "", "fh", "");
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S.SetEvaluator("fa", fa_eval);
    
    // Primary fields
    ep_list.setName("fb");
    // -- field B and its evaluator
    S.Require<CompositeVector,CompositeVectorSpace>("fb", "", "fb");
    fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fb", fb_eval);

    // Primary fields
    ep_list.setName("fc");
    // -- field B and its evaluator
    S.Require<CompositeVector,CompositeVectorSpace>("fc", "", "fc");
    fc_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fc", fc_eval);

    // Primary fields
    ep_list.setName("fh");
    // -- field B and its evaluator
    S.Require<CompositeVector,CompositeVectorSpace>("fh", "", "fh");
    fh_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fh", fh_eval);

    // Primary fields
    ep_list.setName("fe");
    // -- field B and its evaluator
    S.Require<CompositeVector,CompositeVectorSpace>("fe", "", "fe");
    fe_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fe", fe_eval);

    // Setup fields initialize
    S.Setup();
    S.GetW<CompositeVector>("fb", "fb").PutScalar(2.0);
    S.GetRecordW("fb", "fb").set_initialized();
    S.GetW<CompositeVector>("fc", "fc").PutScalar(3.0);
    S.GetRecordW("fc", "fc").set_initialized();
    S.GetW<CompositeVector>("fe", "fe").PutScalar(4.0);
    S.GetRecordW("fe", "fe").set_initialized();
    S.GetW<CompositeVector>("fh", "fh").PutScalar(5.0);
    S.GetRecordW("fh", "fh").set_initialized();
    S.Initialize();
  }

public:
  State S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector,CompositeVectorSpace>> fb_eval, fe_eval, fh_eval, fc_eval;
};

SUITE(DAG) {
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS) {
    // check initialized properly
    CHECK_CLOSE(2.0, (S.Get<CompositeVector>("fb").ViewComponent("cell",false))(0,0), 1e-12);
    CHECK_CLOSE(5.0, (S.Get<CompositeVector>("fh").ViewComponent("cell",false))(0,0), 1e-12);
	
    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    CHECK(changed);
    CHECK_CLOSE(64.0, (S.Get<CompositeVector>("fa").ViewComponent("cell",false))(0,0), 1e-12);
    
    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "main", "fb", "");
    CHECK(changed);
    CHECK_CLOSE(2.0, (S.GetDerivative<CompositeVector>("fa", "", "fb", "").ViewComponent("cell",false))(0,0), 1e-12);
  }
}
