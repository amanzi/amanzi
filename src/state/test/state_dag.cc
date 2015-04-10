/*
  This is the state component of the Amanzi code. 

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
#include "primary_variable_field_evaluator.hh"
#include "secondary_variable_field_evaluator.hh"
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
*/

class AEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  AEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("fa");
    dependencies_.insert(std::string("fb"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {};

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fb_c = *S->GetFieldData("fb")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 2 * fb_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fb_c = *S->GetFieldData("fb")->ViewComponent("cell");

    if (wrt_key == "fb") {
      int ncells = result->size("cell", false);
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 2.0;
      }
    }
  }
};


class make_state {
 public:
  make_state() {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory meshfac(comm);
    mesh = meshfac(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list, ep_list;
    es_list.sublist("VerboseObject").set<std::string>("Verbosity Level", "extreme");

    // Secondary fields
    // -- Field A and its evaluator
    S->RequireField("fa", "fa")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S->SetFieldEvaluator("fa", fa_eval);

    // Primary fields
    // -- field B and its evaluator
    S->RequireField("fb", "fb")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    ep_list.set<std::string>("evaluator name", "fb");
    fb_eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(ep_list));
    S->SetFieldEvaluator("fb", fb_eval);

    // Setup fields and marked as initialized
    S->Setup();
    S->GetField("fa", "fa")->set_initialized();
    S->GetField("fb", "fb")->set_initialized();
    S->Initialize();
  }
  ~make_state() { delete comm; }

 public:
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<State> S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> fb_eval;
};


SUITE(DAG) {
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS) {
    S->GetFieldData("fb", "fb")->PutScalar(2.0);

    fb_eval->SetFieldAsChanged(S.ptr());
    fa_eval->HasFieldChanged(S.ptr(), "main");
    const Epetra_MultiVector& fa = *S->GetFieldData("fa")->ViewComponent("cell");
    CHECK_CLOSE(4.0, fa[0][0], 1e-12);

    fa_eval->HasFieldDerivativeChanged(S.ptr(), "fa", "fb");
    const Epetra_MultiVector& dfa_dfb = *S->GetFieldData("dfa_dfb")->ViewComponent("cell");
    CHECK_CLOSE(2.0, dfa_dfb[0][0], 1e-12);
  }
}

