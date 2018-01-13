/*
  State

  Authors: Ethan Coon

  An improved version of diffusion evaluator forward apply, which builds out
  better evaluators that place lower burden on the PK by better
  EnsureCompatability methods.

  The goal of this is to remove DATA model from PK.
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "TensorVector.hh"
#include "Op_Factory.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Face_Cell.hh"
#include "Operator_Factory.hh"
#include "Operator.hh"
#include "BCs_Factory.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFactory.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "evaluator/EvaluatorSecondary.hh"
#include "evaluator/EvaluatorSecondaries.hh"
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
    double cv = s.GetMesh()->cell_volume(0);
    auto& b = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    b.ViewComponent("cell", false)->PutScalar(-4.*cv);
  }
};

class XIndependent : public EvaluatorIndependent<CompositeVector,CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector,CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new XIndependent(*this)); };

 protected:
  virtual void Update_(State& s) override {
    auto& x = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    auto& x_c = *x.ViewComponent("cell",false);
    for (int c=0; c!=x_c.MyLength(); ++c) {
      AmanziGeometry::Point cc = x.Mesh()->cell_centroid(c);
      x_c[0][c] = cc[0]*cc[0] + cc[1]*cc[1]; // x^2 + y^2      
    }

    if (x.HasComponent("face")) {
      auto& x_f = *x.ViewComponent("face",false);
      for (int f=0; f!=x_f.MyLength(); ++f) {
        AmanziGeometry::Point fc = x.Mesh()->face_centroid(f);
        x_f[0][f] = fc[0]*fc[0] + fc[1]*fc[1]; // x^2 + y^2      
      }
    }
      
    
  }
};

class DiagIndependent : public EvaluatorIndependent<CompositeVector,CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector,CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DiagIndependent(*this)); };

 protected:
  virtual void Update_(State& s) override {
    s.GetW<CompositeVector>(my_key_, my_tag_, my_key_).PutScalar(1.);
  }
};

class KIndependent : public EvaluatorIndependent<TensorVector,TensorVector_Factory> {
 public:
  using EvaluatorIndependent<TensorVector,TensorVector_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new KIndependent(*this)); };

 protected:
  virtual void Update_(State& s) override {
    auto& K = s.GetW<TensorVector>(my_key_, my_tag_, my_key_);
    for (auto& k : K) {
      k.Init(2, 0);
      k.PutScalar(1.0);
    }
  }
};

class BCsIndependent : public EvaluatorIndependent<Operators::BCs,Operators::BCs_Factory> {
 public:
  using EvaluatorIndependent<Operators::BCs,Operators::BCs_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new BCsIndependent(*this)); };

 protected:
  virtual void Update_(State& s) override {
    auto& bcs = s.GetW<Operators::BCs>(my_key_, my_tag_, my_key_);

    auto mesh = bcs.mesh();
    auto kind = bcs.kind();

    auto& model = bcs.bc_model();
    auto& value = bcs.bc_value();
    
    for (auto& bc : model) bc = Operators::OPERATOR_BC_NONE;
    for (auto& val : value) val = 0.;

    // set all exterior faces to dirichlet 0
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells;
    for (int f=0; f!= nfaces_owned; ++f) {
      mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 1) {
        model[f] = Operators::OPERATOR_BC_DIRICHLET;
        auto fc = mesh->face_centroid(f);
        value[f] = fc[0]*fc[0] + fc[1]*fc[1]; // x^2 + y^2
      }
    }    
  }
};


class Evaluator_PDE_DiffusionFV : public EvaluatorSecondaries {
 public:
  Evaluator_PDE_DiffusionFV(Teuchos::ParameterList& plist) :
      EvaluatorSecondaries(plist) {
    tag_ = plist.get<std::string>("tag");

    // my keys
    rhs_key_ = plist.get<std::string>("rhs key");
    local_op_key_ = plist.get<std::string>("local operator key");
    my_keys_.emplace_back(std::make_pair(local_op_key_, tag_));
    my_keys_.emplace_back(std::make_pair(rhs_key_, tag_));

    // dependencies
    tensor_coef_key_ = plist.get<std::string>("tensor coefficient key");
    scalar_coef_key_ = plist.get<std::string>("scalar coefficient key");
    bcs_key_ = plist.get<std::string>("boundary conditions key");
    source_key_ = plist.get<std::string>("source key");
    dependencies_.emplace_back(std::make_pair(tensor_coef_key_, tag_));
    dependencies_.emplace_back(std::make_pair(scalar_coef_key_, tag_));
    dependencies_.emplace_back(std::make_pair(bcs_key_, tag_));
    dependencies_.emplace_back(std::make_pair(source_key_, tag_));
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new Evaluator_PDE_DiffusionFV(*this)); };
  
  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Key& wrt_tag) const override {
    return false;
  }
    
  virtual void EnsureCompatibility(State& S) override {
    // require the rhs
    auto& rhs_fac = S.Require<CompositeVector,CompositeVectorSpace>(rhs_key_, tag_, rhs_key_);
    if (rhs_fac.Mesh().get()) {
      // we have a mesh for the RHS, so we can create a diffusion op to get the schema
      Operators::PDE_DiffusionFactory diff_fac;
      auto diff = diff_fac.Create(plist_, rhs_fac.Mesh());

      // now we can set up the local op
      auto& lop_fac = S.Require<Operators::Op,Operators::Op_Factory>(local_op_key_, tag_, local_op_key_);
      lop_fac.set_mesh(rhs_fac.Mesh());
      Operators::Schema schema(diff->schema_dofs());
      lop_fac.set_schema(schema);

      // push schema to the rhs cvs
      CompositeVectorSpace cvs = Operators::cvsFromSchema(schema, rhs_fac.Mesh());
      rhs_fac.Update(cvs);

      // require scalar coef on the space required by little_k option of operator
      S.Require<CompositeVector,CompositeVectorSpace>(scalar_coef_key_, tag_)
          .Update(diff->little_k_space());

      // require bcs
      auto& bc_fac = S.Require<Operators::BCs,Operators::BCs_Factory>(bcs_key_, tag_);
      bc_fac.set_mesh(rhs_fac.Mesh());
      bc_fac.set_kind(AmanziMesh::FACE);
      bc_fac.set_type(Operators::SCHEMA_DOFS_SCALAR);

      // require source
      auto& src_fac = S.Require<CompositeVector,CompositeVectorSpace>(source_key_, tag_);
      src_fac.SetMesh(rhs_fac.Mesh())
          ->SetGhosted()
          ->AddComponent("cell", AmanziMesh::CELL, 1);

      // require tensors on cells
      auto& K_fac = S.Require<TensorVector,TensorVector_Factory>(tensor_coef_key_, tag_);
      K_fac.set_map(src_fac);
      K_fac.set_rank(0); // REMOVE ME!!!!!!!! GET FROM INDEPENDENT VEC, ONLY FOR TESTING!


    }
  }


  virtual void Update_(State& S) override {
    auto A_rhs = S.GetPtrW<CompositeVector>(rhs_key_, tag_, rhs_key_);
    *A_rhs->ViewComponent("cell", false) = *S.Get<CompositeVector>(source_key_, tag_).ViewComponent("cell", false);
    auto A_lop = S.GetPtrW<Operators::Op>(local_op_key_, tag_, local_op_key_);

    // create the global operator
    Operators::Operator_Factory global_op_fac;
    global_op_fac.set_mesh(A_rhs->Mesh());
    global_op_fac.set_cvs(A_rhs->Map(), A_rhs->Map());

    auto global_op_unique = global_op_fac.Create();
    // need to figure out a way to move unique_ptr into rcp
    // In the mean time, i don't think this will break the world.
    auto global_op = Teuchos::rcpFromRef(*global_op_unique);

    global_op->set_rhs(A_rhs);

    Operators::PDE_DiffusionFactory diff_fac;
    auto pde = diff_fac.Create(plist_, global_op);
    pde->set_local_operator(A_lop);

    Teuchos::RCP<const Operators::BCs> bcs = S.GetPtr<Operators::BCs>(bcs_key_, tag_);
    pde->SetBCs(bcs, bcs);

    const auto& K = S.Get<TensorVector>(tensor_coef_key_, tag_);
    Teuchos::RCP<const std::vector<WhetStone::Tensor>> Kdata = Teuchos::rcpFromRef(K.data);
    pde->SetTensorCoefficient(Kdata);

    // at least this is const!
    Teuchos::RCP<const CompositeVector> kr = S.GetPtr<CompositeVector>(scalar_coef_key_, tag_);
    pde->SetScalarCoefficient(kr, Teuchos::null);

    // compute local ops
    pde->UpdateMatrices(Teuchos::null, Teuchos::null);
    pde->ApplyBCs(true, true);
  }


  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override {
    ASSERT(0);
  }
  
  
 protected:
  Key tag_;
  Key rhs_key_, local_op_key_;
  Key tensor_coef_key_, scalar_coef_key_, source_key_;
  Key bcs_key_;

};


void test(const std::string& discretization) {
    auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory meshfac(comm);
    auto mesh = meshfac(-1.0, -1.0, 1.0, 1.0, 128,128);

    State S;
    S.RegisterDomainMesh(mesh);

    Teuchos::ParameterList es_list;
    es_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    es_list.setName("my_op");

    // NOTE: still need to require evaluators because factory is not tested
    // yet.  Just testing remove of DATA needs.

    // require independent evaluator for x
    Teuchos::ParameterList xe_list;
    xe_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    xe_list.setName("x");
    auto x_eval = Teuchos::rcp(new XIndependent(xe_list));
    S.SetEvaluator("x", x_eval);

    // require independent evaluator for source term b
    Teuchos::ParameterList be_list;
    be_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    be_list.setName("b");
    auto b_eval = Teuchos::rcp(new BIndependent(be_list));
    S.SetEvaluator("b", b_eval);

    // require vector and independent evaluator for Tensor
    Teuchos::ParameterList Ke_list;
    Ke_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    Ke_list.setName("K");
    auto K_eval = Teuchos::rcp(new KIndependent(Ke_list));
    S.SetEvaluator("K", K_eval);
    
    // require vector and independent evaluator for kr (on faces!)
    Teuchos::ParameterList kre_list;
    kre_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    kre_list.setName("k_relative");
    auto kr_eval = Teuchos::rcp(new DiagIndependent(kre_list));
    S.SetEvaluator("k_relative", kr_eval);
    
    // require boundary conditions
    Teuchos::ParameterList bce_list;
    bce_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    bce_list.setName("bcs");
    auto bc_eval = Teuchos::rcp(new BCsIndependent(bce_list));
    S.SetEvaluator("bcs", bc_eval);
    

    // require the local operator and rhs
    Teuchos::ParameterList Ae_list;
    Ae_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    Ae_list.setName("A_local");
    Ae_list.set("tag", "");
    Ae_list.set("local operator key", "A_local");
    Ae_list.set("rhs key", "A_rhs");
    Ae_list.set("tensor coefficient key", "K");
    Ae_list.set("scalar coefficient key", "k_relative");
    Ae_list.set("boundary conditions key", "bcs");
    Ae_list.set("source key", "b");
    Ae_list.set("discretization primary", discretization);
    auto A_eval = Teuchos::rcp(new Evaluator_PDE_DiffusionFV(Ae_list));
    S.SetEvaluator("A_local", A_eval);
    S.SetEvaluator("A_rhs", A_eval);
    
    // require secondary evaluator for r = Ax - b
    Teuchos::ParameterList re_list;
    re_list.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    re_list.set("rhs key", "A_rhs");
    re_list.set("x key", "x");
    re_list.set("local operator keys", Teuchos::Array<std::string>(1,"A_local"));
    re_list.setName("residual");
    auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
    S.SetEvaluator("residual", r_eval);

    // The actual schema has to get set somewhere!  I don't know how it would
    // get set otherwise... lol
    S.Require<CompositeVector,CompositeVectorSpace>("residual", "").SetMesh(mesh)
        ->SetGhosted(true);
    //        ->SetComponent("cell", AmanziMesh::CELL, 1);
    
    // Setup fields and marked as initialized.  Note: USER CODE SHOULD NOT DO IT THIS WAY!
    S.Setup();
    S.Initialize();

    // Update residual
    int updated = S.GetEvaluator("residual")->Update(S, "pk");
    CHECK(updated);

    // b - Ax
    double error(0.);
    auto& r = S.Get<CompositeVector>("residual","");
    r.NormInf(&error);
    std::cout << "Error = " << error << std::endl;
    CHECK(error != 0.0);
    CHECK_CLOSE(0.0, error, 1.e-3);
}
  

SUITE(EVALUATOR_ON_OP) {

  // Apply a non-diagonal operator, including boundary conditions
  TEST(OP_APPLY_DIFFUSION) {
    test("fv: default");
  }


  TEST(OP_APPLY_DIFFUSION_MFD) {
    test("mfd: two-point flux approximation");
  }
  
}
