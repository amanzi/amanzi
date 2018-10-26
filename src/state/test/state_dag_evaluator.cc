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

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "Executor.hh"
#include "state_dag_model.cc"

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
    auto result_c = results[0]->ViewComponent("cell",false);
    auto fb_c = S.Get<CompositeVector>("fb").ViewComponent("cell",false);
    auto fc_c = S.Get<CompositeVector>("fc").ViewComponent("cell",false);
    auto fe_c = S.Get<CompositeVector>("fe").ViewComponent("cell",false);
    auto fh_c = S.Get<CompositeVector>("fh").ViewComponent("cell",false);
    
    AModel<AmanziDefaultDevice> model(result_c, fb_c, fc_c,fe_c,fh_c, plist_);
    
    ExecuteModel("A",model,0,result_c.extent(0));
    ExecuteModel("A",model::dAdB,0,result_c.extent(0));
    ExecuteModel("A",model::dAdC,0,result_c.extent(0));
    ExecuteModel("A",model::dAdE,0,result_c.extent(0));
    ExecuteModel("A",model::dAdH,0,result_c.extent(0));

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
    //    auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    //AmanziMesh::MeshFactory fac(comm);
    MeshFactory meshfac(comm);
    auto mesh = meshfac(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // create a state
    //State S;
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
    //    S.RequireDerivative<CompositeVector,CompositeVectorSpace>("fa", "", "fg", "");
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S.SetEvaluator("fa", fa_eval);
    
    // --  C and its evaluator
    es_list.setName("fc");
    S.Require<CompositeVector,CompositeVectorSpace>("fc", "", "fc")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    fc_eval = Teuchos::rcp(new CEvaluator(es_list));
    S.SetEvaluator("fc", fc_eval);
    /*
    // --  D and its evaluator
    es_list.setName("fd");
    S.Require<CompositeVector,CompositeVectorSpace>("fd", "", "fd")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    fd_eval = Teuchos::rcp(new DEvaluator(es_list));
    S.SetEvaluator("fd", fd_eval);
    */
    // --  E and its evaluator
    es_list.setName("fe");
    S.Require<CompositeVector,CompositeVectorSpace>("fe", "", "fe")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    //    S.RequireDerivative<CompositeVector,CompositeVectorSpace>("fe", "", "fg", "");
    fe_eval = Teuchos::rcp(new EEvaluator(es_list));
    S.SetEvaluator("fe", fe_eval);
    

    /*
    // --  F and its evaluator
    es_list.setName("ff");
    S.Require<CompositeVector,CompositeVectorSpace>("ff", "", "ff")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    ff_eval = Teuchos::rcp(new FEvaluator(es_list));
    S.SetEvaluator("ff", ff_eval);
    */
    // --  H and its evaluator
    es_list.setName("fh");
    S.Require<CompositeVector,CompositeVectorSpace>("fh", "", "fh")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    fh_eval = Teuchos::rcp(new HEvaluator(es_list));
    S.SetEvaluator("fh", fh_eval);

    // Primary fields
    ep_list.setName("fb");
    // -- field B and its evaluator
    S.Require<CompositeVector,CompositeVectorSpace>("fb", "", "fb")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    fb_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fb", fb_eval);
    /*
    // -- field G and its evaluator
    ep_list.setName("fg");
    S.Require<CompositeVector,CompositeVectorSpace>("fg", "", "fg")
      .SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    fg_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(ep_list));
    S.SetEvaluator("fg", fg_eval);
    */
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
    //S.GetW<CompositeVector>("fg", "fg").PutScalar(3.0);
    //S.GetRecordW("fg", "fg").set_initialized();
    S.Initialize();

  }

public:
  State S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<CEvaluator> fc_eval;
  //  Teuchos::RCP<DEvaluator> fd_eval;
  Teuchos::RCP<EEvaluator> fe_eval;
  //  Teuchos::RCP<FEvaluator> ff_eval;
  Teuchos::RCP<HEvaluator> fh_eval;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector,CompositeVectorSpace>> fb_eval, fg_eval;
};

SUITE(DAG) {
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS) {
    // check initialized properly
    CHECK_CLOSE(2.0, (*S.Get<CompositeVector>("fb").ViewComponent("cell",false))[0][0], 1e-12);
    //--CHECK_CLOSE(3.0, (*S.Get<CompositeVector>("fg").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK_CLOSE(5.0, (*S.Get<CompositeVector>("fh").ViewComponent("cell",false))[0][0], 1e-12);
	
    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    CHECK_CLOSE(64.0, (*S.Get<CompositeVector>("fa").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK(changed);

    
    // check intermediate steps got updated too
    //--CHECK_CLOSE(6.0, (*S.Get<CompositeVector>("fd").ViewComponent("cell", false))[0][0], 1e-12);
    
    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fb", "");
    CHECK_CLOSE(2.0, (*S.GetDerivative<CompositeVector>("fa", "", "fb", "").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK(changed);
    /*
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", "");
    CHECK_CLOSE(8640.0, (*S.GetDerivative<CompositeVector>("fa", "", "fg", "").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK(changed);
    
    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    changed = fe_eval->UpdateDerivative(S, "fe", "fg", "");
    CHECK_CLOSE(24.0, (*S.GetDerivative<CompositeVector>("fe", "", "fg", "").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK(changed);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", "");
    CHECK_CLOSE(8640.0, (*S.GetDerivative<CompositeVector>("fa", "", "fg", "").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK(!changed);
    
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetChanged();
    changed = fa_eval->UpdateDerivative(S, "fa", "fg", "");
    CHECK_CLOSE(8640.0, (*S.GetDerivative<CompositeVector>("fa", "", "fg", "").ViewComponent("cell",false))[0][0], 1e-12);
    CHECK(changed);
    */
  }
}
