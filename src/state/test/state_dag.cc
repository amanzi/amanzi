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
class AEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  AEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("fa");
    dependencies_.insert(std::string("fb"));
    dependencies_.insert(std::string("fc"));
    dependencies_.insert(std::string("fe"));
    dependencies_.insert(std::string("fh"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const { return Teuchos::null; }

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fb_c = *S->GetFieldData("fb")->ViewComponent("cell");
    const Epetra_MultiVector& fc_c = *S->GetFieldData("fc")->ViewComponent("cell");
    const Epetra_MultiVector& fe_c = *S->GetFieldData("fe")->ViewComponent("cell");
    const Epetra_MultiVector& fh_c = *S->GetFieldData("fh")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 2 * fb_c[0][c] + fc_c[0][c] * fe_c[0][c] * fh_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fc_c = *S->GetFieldData("fc")->ViewComponent("cell");
    const Epetra_MultiVector& fe_c = *S->GetFieldData("fe")->ViewComponent("cell");
    const Epetra_MultiVector& fh_c = *S->GetFieldData("fh")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    if (wrt_key == "fb") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 2.0;
      }
    } else if (wrt_key == "fc") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = fe_c[0][c] * fh_c[0][c];
      }
    } else if (wrt_key == "fe") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = fc_c[0][c] * fh_c[0][c];
      }
    } else if (wrt_key == "fh") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = fc_c[0][c] * fe_c[0][c];
      }
    }
  }
};


/* ******************************************************************
* Equation C = 2*D + G
****************************************************************** */
class CEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  CEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("fc");
    dependencies_.insert(std::string("fd"));
    dependencies_.insert(std::string("fg"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const { return Teuchos::null; }

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fd_c = *S->GetFieldData("fd")->ViewComponent("cell");
    const Epetra_MultiVector& fg_c = *S->GetFieldData("fg")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 2 * fd_c[0][c] + fg_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    int ncells = result->size("cell", false);

    if (wrt_key == "fd") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 2.0;
      }
    } else if (wrt_key == "fg") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 1.0;
      }
    }
  }
};


/* ******************************************************************
* Equation D = 2*G
****************************************************************** */
class DEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  DEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("fd");
    dependencies_.insert(std::string("fg"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const { return Teuchos::null; }

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fg_c = *S->GetFieldData("fg")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 2 * fg_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    int ncells = result->size("cell", false);

    if (wrt_key == "fg") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 2.0;
      }
    }
  }
};


/* ******************************************************************
* Equation E = D*F
****************************************************************** */
class EEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  EEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("fe");
    dependencies_.insert(std::string("fd"));
    dependencies_.insert(std::string("ff"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const { return Teuchos::null; }

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fd_c = *S->GetFieldData("fd")->ViewComponent("cell");
    const Epetra_MultiVector& ff_c = *S->GetFieldData("ff")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = fd_c[0][c] * ff_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fd_c = *S->GetFieldData("fd")->ViewComponent("cell");
    const Epetra_MultiVector& ff_c = *S->GetFieldData("ff")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    if (wrt_key == "fd") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = ff_c[0][c];
      }
    } else if (wrt_key == "ff") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = fd_c[0][c];
      }
    }
  }
};


/* ******************************************************************
* Equation F = 2*G
****************************************************************** */
class FEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  FEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("ff");
    dependencies_.insert(std::string("fg"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const { return Teuchos::null; }

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& fg_c = *S->GetFieldData("fg")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 2 * fg_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    int ncells = result->size("cell", false);

    if (wrt_key == "fg") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 2.0;
      }
    }
  }
};


/* ******************************************************************
* Equation H = 2*F
****************************************************************** */
class HEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  HEvaluator(Teuchos::ParameterList& plist) : 
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = std::string("fh");
    dependencies_.insert(std::string("ff"));
  }

  virtual Teuchos::RCP<FieldEvaluator> Clone() const { return Teuchos::null; }

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    const Epetra_MultiVector& ff_c = *S->GetFieldData("ff")->ViewComponent("cell");

    int ncells = result->size("cell", false);
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 2 * ff_c[0][c];
    }
  }

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");
    int ncells = result->size("cell", false);

    if (wrt_key == "ff") {
      for (int c = 0; c != ncells; ++c) {
        result_c[0][c] = 2.0;
      }
    }
  }
};


class make_state {
 public:
  make_state() {
    comm = Amanzi::getDefaultComm();
    MeshFactory meshfactory(comm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    S = Teuchos::rcp(new State());
    S->RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list, ep_list;
    es_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");

    // Secondary fields
    // -- Field A and its evaluator
    S->RequireField("fa", "fa")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    fa_eval = Teuchos::rcp(new AEvaluator(es_list));
    S->SetFieldEvaluator("fa", fa_eval);

    // -- Field C and its evaluator
    S->RequireField("fc", "fc")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    fc_eval = Teuchos::rcp(new CEvaluator(es_list));
    S->SetFieldEvaluator("fc", fc_eval);

    // -- Field D and its evaluator
    S->RequireField("fd", "fd")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    fd_eval = Teuchos::rcp(new DEvaluator(es_list));
    S->SetFieldEvaluator("fd", fd_eval);

    // -- Field E and its evaluator
    S->RequireField("fe", "fe")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    fe_eval = Teuchos::rcp(new EEvaluator(es_list));
    S->SetFieldEvaluator("fe", fe_eval);

    // -- Field F and its evaluator
    S->RequireField("ff", "ff")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    ff_eval = Teuchos::rcp(new FEvaluator(es_list));
    S->SetFieldEvaluator("ff", ff_eval);

    // -- Field H and its evaluator
    S->RequireField("fh", "fh")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    fh_eval = Teuchos::rcp(new HEvaluator(es_list));
    S->SetFieldEvaluator("fh", fh_eval);

    // Primary fields
    // -- field B and its evaluator
    S->RequireField("fb", "fb")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    ep_list.set<std::string>("evaluator name", "fb");
    fb_eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(ep_list));
    S->SetFieldEvaluator("fb", fb_eval);

    // -- field G and its evaluator
    S->RequireField("fg", "fg")->SetMesh(mesh)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    ep_list.set<std::string>("evaluator name", "fg");
    fg_eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(ep_list));
    S->SetFieldEvaluator("fg", fg_eval);

    // Setup fields and marked as initialized
    S->Setup();
    S->GetField("fb", "fb")->set_initialized();
    S->GetField("fg", "fg")->set_initialized();
    S->Initialize();
  }
  ~make_state() { }

 public:
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<State> S;
  Teuchos::RCP<AEvaluator> fa_eval;
  Teuchos::RCP<CEvaluator> fc_eval;
  Teuchos::RCP<DEvaluator> fd_eval;
  Teuchos::RCP<EEvaluator> fe_eval;
  Teuchos::RCP<FEvaluator> ff_eval;
  Teuchos::RCP<HEvaluator> fh_eval;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> fb_eval, fg_eval;
};


SUITE(DAG) {
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS) {
    // set primary fields
    S->GetFieldData("fb", "fb")->PutScalar(2.0);
    fb_eval->SetFieldAsChanged(S.ptr());

    S->GetFieldData("fg", "fg")->PutScalar(3.0);
    fg_eval->SetFieldAsChanged(S.ptr());

    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    fa_eval->HasFieldChanged(S.ptr(), "main");
    const Epetra_MultiVector& fa = *S->GetFieldData("fa")->ViewComponent("cell");
    CHECK_CLOSE(6484.0, fa[0][0], 1e-12);

    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    fa_eval->HasFieldDerivativeChanged(S.ptr(), "fa", "fb");
    const Epetra_MultiVector& dfa_dfb = *S->GetFieldData(Keys::getDerivKey("fa", "fb"))->ViewComponent("cell");
    CHECK_CLOSE(2.0, dfa_dfb[0][0], 1e-12);

    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fa_eval->HasFieldDerivativeChanged(S.ptr(), "fa", "fg");
    const Epetra_MultiVector& dfa_dfg = *S->GetFieldData(Keys::getDerivKey("fa","fg"))->ViewComponent("cell");
    CHECK_CLOSE(8640.0, dfa_dfg[0][0], 1e-12);

    // calculate dE/dD is not well defined.
    // We keep the code which actually work correctly.
    /*
    fg_eval->SetFieldAsChanged(S.ptr());
    std::cout << "Calculate derivative of field E wrt field D:" << std::endl;
    fe_eval->HasFieldDerivativeChanged(S.ptr(), "fe", "fd");
    const Epetra_MultiVector& dfe_dfd = *S->GetFieldData("dfe_dfd")->ViewComponent("cell");
    CHECK_CLOSE(6.0, dfe_dfd[0][0], 1e-12);
    */

    // calculate dE/dB: This is really strange
    /*
    std::cout << "Calculate derivative of field E wrt field B:" << std::endl;
    fe_eval->HasFieldDerivativeChanged(S.ptr(), "fe", "fb");
    const Epetra_MultiVector& dfe_dfb = *S->GetFieldData("dfe_dfb")->ViewComponent("cell");
    CHECK_CLOSE(0.0, dfe_dfb[0][0], 1e-12);
    */

    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    fe_eval->HasFieldDerivativeChanged(S.ptr(), "fe", "fg");
    const Epetra_MultiVector& dfe_dfg = *S->GetFieldData(Keys::getDerivKey("fe","fg"))->ViewComponent("cell");
    CHECK_CLOSE(24.0, dfe_dfg[0][0], 1e-12);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetFieldAsChanged(S.ptr());
    fa_eval->HasFieldDerivativeChanged(S.ptr(), "fa", "fg");
    const Epetra_MultiVector& dfa_dfg2 = *S->GetFieldData(Keys::getDerivKey("fa","fg"))->ViewComponent("cell");
    CHECK_CLOSE(8640.0, dfa_dfg2[0][0], 1e-12);
  }
}

