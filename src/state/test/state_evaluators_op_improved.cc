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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "BCs_Factory.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Face_Cell.hh"
#include "Op_Factory.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFactory.hh"
#include "State.hh"
#include "TensorVector.hh"
#include "evaluator/EvaluatorAlgebraic.hh"
#include "evaluator/EvaluatorIndependent.hh"
#include "evaluator/EvaluatorPrimary.hh"
#include "evaluator/EvaluatorSecondaries.hh"
#include "evaluator/EvaluatorSecondary.hh"
#include "evaluator/Evaluator_OperatorApply.hh"
#include "evaluator/Evaluator_PDE_Diffusion.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

class BIndependent
    : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
public:
  using EvaluatorIndependent<CompositeVector,
                             CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new BIndependent(*this));
  };

protected:
  virtual void Update_(State &s) override {
    double cv = s.GetMesh()->cell_volume(0);
    auto &b = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    b.ViewComponent("cell", false)->PutScalar(-4. * cv);
  }
};

class XIndependent
    : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
public:
  using EvaluatorIndependent<CompositeVector,
                             CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new XIndependent(*this));
  };

protected:
  virtual void Update_(State &s) override {
    auto &x = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    auto &x_c = *x.ViewComponent("cell", false);
    for (int c = 0; c != x_c.MyLength(); ++c) {
      AmanziGeometry::Point cc = x.Mesh()->cell_centroid(c);
      x_c[0][c] = cc[0] * cc[0] + cc[1] * cc[1]; // x^2 + y^2
    }

    if (x.HasComponent("face")) {
      auto &x_f = *x.ViewComponent("face", false);
      for (int f = 0; f != x_f.MyLength(); ++f) {
        AmanziGeometry::Point fc = x.Mesh()->face_centroid(f);
        x_f[0][f] = fc[0] * fc[0] + fc[1] * fc[1]; // x^2 + y^2
      }
    }
  }
};

class DiagIndependent
    : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
public:
  using EvaluatorIndependent<CompositeVector,
                             CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DiagIndependent(*this));
  };

protected:
  virtual void Update_(State &s) override {
    s.GetW<CompositeVector>(my_key_, my_tag_, my_key_).PutScalar(1.);
  }
};

class KIndependent
    : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
public:
  using EvaluatorIndependent<TensorVector,
                             TensorVector_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new KIndependent(*this));
  };

protected:
  virtual void Update_(State &s) override {
    auto &K = s.GetW<TensorVector>(my_key_, my_tag_, my_key_);
    for (auto &k : K) {
      // note hard coded rank here, the actual tensor independent eval should
      // set this in consistent, and it should get init'ed somewhere during that
      // process, not now!
      k.Init(2, 1);
      k.PutScalar(1.0);
    }
  }
};

class BCsIndependent
    : public EvaluatorIndependent<Operators::BCs, Operators::BCs_Factory> {
public:
  using EvaluatorIndependent<Operators::BCs,
                             Operators::BCs_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new BCsIndependent(*this));
  };

protected:
  virtual void Update_(State &s) override {
    auto &bcs = s.GetW<Operators::BCs>(my_key_, my_tag_, my_key_);

    auto mesh = bcs.mesh();
    auto kind = bcs.kind();

    auto &model = bcs.bc_model();
    auto &value = bcs.bc_value();

    for (auto &bc : model)
      bc = Operators::OPERATOR_BC_NONE;
    for (auto &val : value)
      val = 0.;

    // set all exterior faces to dirichlet 0
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells;
    for (int f = 0; f != nfaces_owned; ++f) {
      mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 1) {
        model[f] = Operators::OPERATOR_BC_DIRICHLET;
        auto fc = mesh->face_centroid(f);
        value[f] = fc[0] * fc[0] + fc[1] * fc[1]; // x^2 + y^2
      }
    }
  }
};

class AddAlgebraic
    : public EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace> {
public:
  AddAlgebraic(Teuchos::ParameterList &plist)
      : EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace>(plist),
        coefs_(plist.get<Teuchos::Array<double>>("coefficients")) {
    if (coefs_.size() != dependencies_.size()) {
      ASSERT(0);
    }
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new AddAlgebraic(*this));
  };

  virtual void Evaluate_(const State &S, CompositeVector &result) override {
    int i = 0;
    for (const auto &dep : dependencies_) {
      result.Update(coefs_[i], S.Get<CompositeVector>(dep.first, dep.second),
                    1.);
      ++i;
    }
  }
  virtual void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                          const Key &wrt_tag,
                                          CompositeVector &result) override {
    int i = 0;
    for (const auto &dep : dependencies_) {
      if (dep.first == wrt_key && dep.second == wrt_tag) {
        const auto &x = S.Get<CompositeVector>(dep.first, dep.second);
        for (const auto &comp : x)
          result.ViewComponent(comp, false)->PutScalar(coefs_[i]);
        return;
      }
      ++i;
    }
  }

protected:
  Teuchos::Array<double> coefs_;
};

void test(const std::string &discretization) {
  auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  MeshFactory meshfac(comm);
  auto mesh = meshfac(-1.0, -1.0, 1.0, 1.0, 128, 128);

  State S;
  S.RegisterDomainMesh(mesh);

  Teuchos::ParameterList es_list;
  es_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  es_list.setName("my_op");

  // NOTE: still need to require evaluators because factory is not tested
  // yet.  Just testing remove of DATA needs.

  // require independent evaluator for x
  Teuchos::ParameterList xe_list;
  xe_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  xe_list.setName("x");
  auto x_eval = Teuchos::rcp(new XIndependent(xe_list));
  S.SetEvaluator("x", x_eval);

  // require independent evaluator for source term b
  S.Require<CompositeVector, CompositeVectorSpace>("b", "")
      .SetMesh(mesh)
      ->SetGhosted(false)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList be_list;
  be_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  be_list.setName("b");
  auto b_eval = Teuchos::rcp(new BIndependent(be_list));
  S.SetEvaluator("b", b_eval);

  // require vector and independent evaluator for Tensor
  Teuchos::ParameterList Ke_list;
  Ke_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  Ke_list.setName("K");
  auto K_eval = Teuchos::rcp(new KIndependent(Ke_list));
  S.SetEvaluator("K", K_eval);

  // require vector and independent evaluator for kr (on faces!)
  Teuchos::ParameterList kre_list;
  kre_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  kre_list.setName("k_relative");
  auto kr_eval = Teuchos::rcp(new DiagIndependent(kre_list));
  S.SetEvaluator("k_relative", kr_eval);

  // require boundary conditions
  Teuchos::ParameterList bce_list;
  bce_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  bce_list.setName("bcs");
  auto bc_eval = Teuchos::rcp(new BCsIndependent(bce_list));
  S.SetEvaluator("bcs", bc_eval);

  // require the local operator and rhs
  Teuchos::ParameterList Ae_list;
  Ae_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  Ae_list.setName("A_local");
  Ae_list.set("tag", "");
  Ae_list.set("rhs key", "A_rhs");
  Ae_list.set("local operator key", "A_local");
  Ae_list.set("tensor coefficient key", "K");
  Ae_list.set("scalar coefficient key", "k_relative");
  Ae_list.set("boundary conditions key", "bcs");
  Ae_list.set("operator argument key", "x");
  Ae_list.set("discretization primary", discretization);
  auto A_eval = Teuchos::rcp(new Evaluator_PDE_Diffusion(Ae_list));
  S.SetEvaluator("A_local", A_eval);
  S.SetEvaluator("A_rhs", A_eval);

  // require secondary evaluator for r = Ax - b
  Teuchos::ParameterList re_list;
  re_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  re_list.set("rhs key", "A_rhs");
  re_list.set("x key", "x");
  re_list.set("local operator keys", Teuchos::Array<std::string>(1, "A_local"));
  re_list.setName("rhs_m_Ax");
  auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
  S.SetEvaluator("rhs_m_Ax", r_eval);
  S.Require<CompositeVector, CompositeVectorSpace>("rhs_m_Ax", "")
      .SetMesh(mesh)
      ->SetGhosted(true);

  // residual = rhs-Ax + b
  Teuchos::ParameterList res_list("residual");
  res_list.sublist("verbose object")
      .set<std::string>("verbosity level", "extreme");
  res_list.set("tag", "");
  res_list.set("dependencies", Teuchos::Array<std::string>(
                                   std::vector<std::string>{"rhs_m_Ax", "b"}));
  res_list.set("dependency tags",
               Teuchos::Array<std::string>(std::vector<std::string>{"", ""}));
  res_list.set("coefficients",
               Teuchos::Array<double>(std::vector<double>(2, 1.0)));
  res_list.set("consistency policy", "take from child: union");
  S.Require<CompositeVector, CompositeVectorSpace>("residual", "")
      .SetMesh(mesh)
      ->SetGhosted(true);
  auto res_eval = Teuchos::rcp(new AddAlgebraic(res_list));
  S.SetEvaluator("residual", res_eval);

  // Setup fields and marked as initialized.  Note: USER CODE SHOULD NOT DO IT
  // THIS WAY!
  S.Setup();
  S.Initialize();

  // Update residual
  int updated = S.GetEvaluator("residual")->Update(S, "pk");
  CHECK(updated);

  // b - Ax
  double error(0.);
  auto &r = S.Get<CompositeVector>("residual", "");
  r.NormInf(&error);
  std::cout << "Error = " << error << std::endl;
  CHECK(error != 0.0);
  CHECK_CLOSE(0.0, error, 1.e-3);
}

SUITE(EVALUATOR_ON_OP) {

  // Apply a non-diagonal operator, including boundary conditions
  TEST(OP_APPLY_DIFFUSION_FV) { test("fv: default"); }

  TEST(OP_APPLY_DIFFUSION_MFD) { test("mfd: two-point flux approximation"); }

  TEST(OP_APPLY_DIFFUSION_NLFV) { test("nlfv: default"); }
}
