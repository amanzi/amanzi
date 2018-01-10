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
#include "Op_Factory.hh"
#include "Op_Cell_Cell.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "evaluator/EvaluatorSecondary.hh"
#include "evaluator/EvaluatorIndependent.hh"
#include "evaluator/Evaluator_OperatorApply.hh"
#include "State.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

class BIndependent : public EvaluatorIndependent<CompositeVector,CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector,CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new BIndependent(*this)); };

 protected:
  virtual void Update_(State& s) override {
    s.GetW<CompositeVector>(my_key_, my_tag_, my_key_).PutScalar(3.);
  }
};


class DiagIndependent : public EvaluatorIndependent<CompositeVector,CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector,CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DiagIndependent(*this)); };

 protected:
  virtual void Update_(State& s) override {
    s.GetW<CompositeVector>(my_key_, my_tag_, my_key_).PutScalar(6.);
  }
};


class Evaluator_PDE_Diagonal : public EvaluatorSecondary<Operators::Op,Operators::Op_Factory<Operators::Op_Cell_Cell>> {
 public:
  using EvaluatorSecondary<Operators::Op,Operators::Op_Factory<Operators::Op_Cell_Cell>>::EvaluatorSecondary;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new Evaluator_PDE_Diagonal(*this)); };

 protected:
  virtual void Evaluate_(const State& s, Operators::Op& op) override {
    *op.diag = *s.Get<CompositeVector>(dependencies_.begin()->first, dependencies_.begin()->second)
                .ViewComponent("cell", false);
  }

  virtual void EvaluatePartialDerivative_(const State& s, const Key& wrt_key, const Key& wrt_tag,
          Operators::Op& op) override {
    ASSERT(0);
  }
  
  
};



SUITE(EVALUATOR_ON_OP) {
  TEST(PRIMARY) {
    auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory meshfac(comm);
    auto mesh = meshfac(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    es_list.setName("my_op");

    // require some op data
    auto& f = S.Require<Operators::Op, Operators::Op_Factory<Operators::Op_Cell_Cell>>("my_op", "", "my_op");
    f.set_mesh(mesh);
    f.set_name("cell");

    // make a primary evaluator for it
    auto op_eval = Teuchos::rcp(new EvaluatorPrimary<Operators::Op, Operators::Op_Factory<Operators::Op_Cell_Cell>>(es_list));
    S.SetEvaluator("my_op", op_eval);

    // Setup fields and marked as initialized.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();
    S.GetW<Operators::Op>("my_op", "", "my_op").diag->PutScalar(3.14);
    S.GetRecordW("my_op", "my_op").set_initialized();
    S.Initialize();
  }


  TEST(OP_APPLY) {
    auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory meshfac(comm);
    auto mesh = meshfac(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    State S;
    S.RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    es_list.setName("my_op");

    // require vector and primary evaluator for x
    S.Require<CompositeVector,CompositeVectorSpace>("x", "").SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList xe_list;
    xe_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    xe_list.setName("x");
    auto x_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(xe_list));
    S.SetEvaluator("x", x_eval);

    // require vector and independent evaluator for b
    S.Require<CompositeVector,CompositeVectorSpace>("b", "").SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList be_list;
    be_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    be_list.setName("b");
    auto b_eval = Teuchos::rcp(new BIndependent(be_list));
    S.SetEvaluator("b", b_eval);
    
    // require vector and independent evaluator for Diag(A)
    S.Require<CompositeVector,CompositeVectorSpace>("Diag", "").SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList de_list;
    de_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    de_list.setName("Diag");
    auto D_eval = Teuchos::rcp(new DiagIndependent(de_list));
    S.SetEvaluator("Diag", D_eval);
    
    // require the local operator and evaluator
    auto& fac = S.Require<Operators::Op,Operators::Op_Factory<Operators::Op_Cell_Cell>>("Alocal", "");
    fac.set_mesh(mesh);
    Teuchos::ParameterList Ae_list;
    Ae_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    Ae_list.setName("Alocal");
    Ae_list.set("dependencies", Teuchos::Array<std::string>(1,"Diag"));
    Ae_list.set("dependency tags are my tag", true);
    auto A_eval = Teuchos::rcp(new Evaluator_PDE_Diagonal(Ae_list));
    S.SetEvaluator("Alocal", A_eval);
    
    // require vector and secondary evaluator for r = Ax - b
    S.Require<CompositeVector,CompositeVectorSpace>("residual", "").SetMesh(mesh)
        ->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList re_list;
    re_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    re_list.set("rhs key", "b");
    re_list.set("x key", "x");
    re_list.set("local operator keys", Teuchos::Array<std::string>(1,"Alocal"));
    re_list.setName("residual");
    auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
    S.SetEvaluator("residual", r_eval);

    // Setup fields and marked as initialized.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();
    S.GetW<CompositeVector>("x", "", "x").PutScalar(1.);
    S.GetRecordW("x", "", "x").set_initialized();
    S.Initialize();

    // Update residual
    CHECK(S.GetEvaluator("residual")->Update(S, "pk"));

    // b - Ax
    CHECK_CLOSE(-3.0, (*S.Get<CompositeVector>("residual","").ViewComponent("cell", false))[0][0], 1.e-10);
    
  }
  

}

