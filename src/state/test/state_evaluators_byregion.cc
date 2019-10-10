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
#include "GeometricModel.hh"
#include "RegionBox.hh"
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
  AModel(double alpha) :
      alpha_(alpha) {}

  void SetViews(VectorView_type<DeviceType> A, cVectorView_type<DeviceType> B) {
    A_ = A;
    B_ = B;
  }
  
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    A_(i) = alpha_ * B_(i);
  }

  class dAdB {};
  KOKKOS_INLINE_FUNCTION void operator() (dAdB, const int i) const {
    A_(i) = alpha_;
  }

private:
  VectorView_type<DeviceType> A_;
  cVectorView_type<DeviceType> B_;
  double alpha_;
};



/* ******************************************************************
 * Equation A = alpha * B
 ****************************************************************** */
// Device Type : GPUs or CPUs??
class AEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
public:
  AEvaluator(Teuchos::ParameterList &plist)
    : EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace>(plist) {
    dependencies_.emplace_back(std::make_pair(Key("fb"), Key()));

    models_.emplace_back("region1", AModel<AmanziDefaultDevice>(2.0));
    models_.emplace_back("region2", AModel<AmanziDefaultDevice>(4.0));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AEvaluator(*this));
  }

  virtual void Evaluate_(const State &S, const std::vector<CompositeVector*> &results) override {
    const AmanziMesh::Mesh& mesh = *results[0]->Mesh();
    auto result_c = results[0]->ViewComponent("cell",0,false);
    auto fb_c = S.Get<CompositeVector>("fb").ViewComponent("cell",0,false);

    for (auto& model : models_) {
      model.second.SetViews(result_c, fb_c);
      Entity_ID_View region_entities;
      mesh.get_set_entities(model.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED, region_entities);
      ExecuteModel("A", model.second, region_entities);
    }
  }

  virtual void EvaluatePartialDerivative_(const State &S, 
					  const Key &wrt_key,
                                          const Key &wrt_tag,
                                          const std::vector<CompositeVector*> &results) override {
    const AmanziMesh::Mesh& mesh = *results[0]->Mesh();
    auto result = results[0]->ViewComponent("cell",0,false);
    auto fb = S.Get<CompositeVector>("fb").ViewComponent("cell",0,false);

    for (auto& model : models_) {
      model.second.SetViews(result, fb);
      Entity_ID_View region_entities;
      mesh.get_set_entities(model.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED, region_entities);
      ExecuteModel("A", model.second, AModel<AmanziDefaultDevice>::dAdB(), region_entities);
    }
  }

 protected:
  std::vector<std::pair<std::string, AModel<AmanziDefaultDevice> > > models_;
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
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3));
    Teuchos::RCP<AmanziGeometry::Region> r1 = Teuchos::rcp(new AmanziGeometry::RegionBox("region1", 0,
            AmanziGeometry::Point(0.0, 0.0, 0.0), AmanziGeometry::Point(2.0,4.0,4.0)));
    Teuchos::RCP<AmanziGeometry::Region> r2 = Teuchos::rcp(new AmanziGeometry::RegionBox("region2", 0,
            AmanziGeometry::Point(2.0, 0.0, 0.0), AmanziGeometry::Point(4.0,4.0,4.0)));
    gm->AddRegion(r1);
    gm->AddRegion(r2);
    
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm, gm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

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
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S.SetEvaluator("fa", fa_eval);
    
    // Primary fields
    ep_list.setName("fb");
    // -- field B and its evaluator
    S.Require<CompositeVector,CompositeVectorSpace>("fb", "", "fb");
    fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fb", fb_eval);

    // Setup fields initialize
    S.Setup();
    S.GetW<CompositeVector>("fb", "fb").putScalar(2.0);
    S.GetRecordW("fb", "fb").set_initialized();
    S.Initialize();
  }

public:
  State S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector,CompositeVectorSpace> > fb_eval;
};

SUITE(DAG) {
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS) {
    // check initialized properly
    CHECK_CLOSE(2.0, (S.Get<CompositeVector>("fb").ViewComponent<AmanziDefaultHost>("cell",false))(0,0), 1e-12);
    CHECK_CLOSE(2.0, (S.Get<CompositeVector>("fb").ViewComponent<AmanziDefaultHost>("cell",false))(7,0), 1e-12);
	
    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    CHECK(changed);
    CHECK_CLOSE(4.0, (S.Get<CompositeVector>("fa").ViewComponent<AmanziDefaultHost>("cell",false))(0,0), 1e-12);
    CHECK_CLOSE(4.0, (S.Get<CompositeVector>("fa").ViewComponent<AmanziDefaultHost>("cell",false))(1,0), 1e-12);
    CHECK_CLOSE(8.0, (S.Get<CompositeVector>("fa").ViewComponent<AmanziDefaultHost>("cell",false))(6,0), 1e-12);
    CHECK_CLOSE(8.0, (S.Get<CompositeVector>("fa").ViewComponent<AmanziDefaultHost>("cell",false))(7,0), 1e-12);
    
    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "main", "fb", "");
    CHECK(changed);

    CHECK_CLOSE(2.0, (S.GetDerivative<CompositeVector>("fa", "", "fb", "").ViewComponent<AmanziDefaultHost>("cell",false))(0,0), 1e-12);
    CHECK_CLOSE(2.0, (S.GetDerivative<CompositeVector>("fa", "", "fb", "").ViewComponent<AmanziDefaultHost>("cell",false))(1,0), 1e-12);
    CHECK_CLOSE(4.0, (S.GetDerivative<CompositeVector>("fa", "", "fb", "").ViewComponent<AmanziDefaultHost>("cell",false))(6,0), 1e-12);
    CHECK_CLOSE(4.0, (S.GetDerivative<CompositeVector>("fa", "", "fb", "").ViewComponent<AmanziDefaultHost>("cell",false))(7,0), 1e-12);

  }
}

