/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

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

// Amanzi
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

// Amanzi::State
#include "EvaluatorIndependent.hh"
#include "Evaluator_OperatorApply.hh"
#include "EvaluatorPrimary.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "EvaluatorSecondary.hh"
#include "Evaluator_PDE_Diffusion.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

class BIndependent : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector, CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new BIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    double cv = s.GetMesh()->getCellVolume(0);
    auto& b = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    b.ViewComponent("cell", false)->PutScalar(-4. * cv);
  }
};

class XIndependent : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector, CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new XIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    auto& x = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    auto& x_c = *x.ViewComponent("cell");
    for (int c = 0; c != x_c.MyLength(); ++c) {
      AmanziGeometry::Point cc = x.Mesh()->getCellCentroid(c);
      x_c[0][c] = cc[0] * cc[0] + cc[1] * cc[1]; // x^2 + y^2
    }

    if (x.HasComponent("face")) {
      auto& x_f = *x.ViewComponent("face");
      for (int f = 0; f != x_f.MyLength(); ++f) {
        AmanziGeometry::Point fc = x.Mesh()->getFaceCentroid(f);
        x_f[0][f] = fc[0] * fc[0] + fc[1] * fc[1]; // x^2 + y^2
      }
    }
  }
};

class DiagIndependent : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector, CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new DiagIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    s.GetW<CompositeVector>(my_key_, my_tag_, my_key_).PutScalar(1.);
  }
};

class KIndependent : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  using EvaluatorIndependent<TensorVector, TensorVector_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new KIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    auto& K = s.GetW<TensorVector>(my_key_, my_tag_, my_key_);
    for (auto& k : K) {
      // note hard coded rank here, the actual tensor independent eval should
      // set this in consistent, and it should get init'ed somewhere during that
      // process, not now!
      k.Init(2, 1);
      k.PutScalar(1.0);
    }
  }
};

class BCsIndependent : public EvaluatorIndependent<Operators::BCs, Operators::BCs_Factory> {
 public:
  using EvaluatorIndependent<Operators::BCs, Operators::BCs_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new BCsIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    auto& bcs = s.GetW<Operators::BCs>(my_key_, my_tag_, my_key_);

    auto mesh = bcs.mesh();

    auto& model = bcs.bc_model();
    auto& value = bcs.bc_value();

    for (auto& bc : model) bc = Operators::OPERATOR_BC_NONE;
    for (auto& val : value) val = 0.;

    // set all exterior faces to dirichlet 0
    int nfaces_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    AmanziMesh::Entity_ID_View cells;
    for (int f = 0; f != nfaces_owned; ++f) {
      cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      if (cells.size() == 1) {
        model[f] = Operators::OPERATOR_BC_DIRICHLET;
        auto fc = mesh->getFaceCentroid(f);
        value[f] = fc[0] * fc[0] + fc[1] * fc[1]; // x^2 + y^2
      }
    }
  }
};


void
test(const std::string& discretization)
{
  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfac(comm);
  auto mesh = meshfac.create(-1.0, -1.0, 1.0, 1.0, 128, 128);

  State S;
  S.RegisterDomainMesh(mesh);

  std::cout << "\n\n========= NEXT test for \"" << discretization << "\"\n";
  Teuchos::ParameterList es_list;
  es_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  es_list.setName("my_op");

  // NOTE: still need to require evaluators because factory is not tested
  // yet.  Just testing remove of DATA needs.

  // require independent evaluator for x
  Teuchos::ParameterList xe_list;
  xe_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  xe_list.setName("x");
  auto x_eval = Teuchos::rcp(new XIndependent(xe_list));
  S.SetEvaluator("x", Tags::DEFAULT, x_eval);

  // require independent evaluator for source term b
  S.Require<CompositeVector, CompositeVectorSpace>("b", Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::ParameterList be_list;
  be_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  be_list.setName("b");
  auto b_eval = Teuchos::rcp(new BIndependent(be_list));
  S.SetEvaluator("b", Tags::DEFAULT, b_eval);

  // require vector and independent evaluator for Tensor
  Teuchos::ParameterList Ke_list;
  Ke_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  Ke_list.setName("K");
  auto K_eval = Teuchos::rcp(new KIndependent(Ke_list));
  S.SetEvaluator("K", Tags::DEFAULT, K_eval);

  // require vector and independent evaluator for kr (on faces!)
  Teuchos::ParameterList kre_list;
  kre_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  kre_list.setName("k_relative");
  auto kr_eval = Teuchos::rcp(new DiagIndependent(kre_list));
  S.SetEvaluator("k_relative", Tags::DEFAULT, kr_eval);

  // require boundary conditions
  Teuchos::ParameterList bce_list;
  bce_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  bce_list.setName("bcs");
  auto bc_eval = Teuchos::rcp(new BCsIndependent(bce_list));
  S.SetEvaluator("bcs", Tags::DEFAULT, bc_eval);

  // require the local operator and rhs
  Teuchos::ParameterList Ae_list;
  Ae_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  Ae_list.setName("A_local")
    .set("tag", "")
    .set("rhs key", "A_rhs")
    .set("local operator key", "A_local")
    .set("tensor coefficient key", "K")
    .set("scalar coefficient key", "k_relative")
    .set("boundary conditions key", "bcs")
    .set("operator argument key", "x")
    .set("discretization primary", discretization);

  auto A_eval = Teuchos::rcp(new Evaluator_PDE_Diffusion(Ae_list));
  S.SetEvaluator("A_local", Tags::DEFAULT, A_eval);
  S.SetEvaluator("A_rhs", Tags::DEFAULT, A_eval);

  // require secondary evaluator for r = Ax - b
  Teuchos::ParameterList re_list("residual");
  re_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  re_list.set("diagonal primary x key", "x");
  re_list.set("diagonal local operators keys", Teuchos::Array<std::string>(1, "A_local"));
  re_list.set("diagonal local operator rhss keys", Teuchos::Array<std::string>(1, "A_rhs"));
  re_list.set("additional rhss keys", Teuchos::Array<std::string>(1, "b"));
  re_list.set("rhs coefficients", Teuchos::Array<double>(1, -1.0));
  re_list.set("tag", "");
  auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
  S.SetEvaluator("residual", Tags::DEFAULT, r_eval);
  S.Require<CompositeVector, CompositeVectorSpace>("residual", Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetGhosted(true);

  // run
  // -- setup and init
  S.Setup();
  S.Initialize();

  // -- update residual
  int updated = S.GetEvaluator("residual").Update(S, "pk");
  CHECK(updated);

  // -- check error
  double error;
  const auto& r = S.Get<CompositeVector>("residual", Tags::DEFAULT);
  r.NormInf(&error);
  std::cout << "Error=" << error << std::endl;
  CHECK_CLOSE(0.0, error, 1.e-3);
}


void
test_inverse(const std::string& discretization)
{
  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfac(comm);
  // auto mesh = meshfac(-1.0, -1.0, 1.0, 1.0, 128, 128);
  auto mesh = meshfac.create(-1.0, -1.0, 1.0, 1.0, 4, 4);

  State S;
  S.RegisterDomainMesh(mesh);

  std::cout << "\n\n========= NEXT test for \"" << discretization << "\"\n";
  Teuchos::ParameterList es_list;
  es_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  es_list.setName("my_op");

  // NOTE: still need to require evaluators because factory is not tested
  // yet.  Just testing remove of DATA needs.

  // require primary evaluator for x
  S.Require<CompositeVector, CompositeVectorSpace>("x", Tags::DEFAULT, "x")
    .SetMesh(mesh)
    ->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::ParameterList xe_list;
  xe_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  xe_list.setName("x");
  auto x_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(xe_list));
  S.SetEvaluator("x", Tags::DEFAULT, x_eval);

  // require independent evaluator for source term b
  S.Require<CompositeVector, CompositeVectorSpace>("b", Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetGhosted(false)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  Teuchos::ParameterList be_list;
  be_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  be_list.setName("b");
  auto b_eval = Teuchos::rcp(new BIndependent(be_list));
  S.SetEvaluator("b", Tags::DEFAULT, b_eval);

  // require vector and independent evaluator for Tensor
  Teuchos::ParameterList Ke_list;
  Ke_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  Ke_list.setName("K");
  auto K_eval = Teuchos::rcp(new KIndependent(Ke_list));
  S.SetEvaluator("K", Tags::DEFAULT, K_eval);

  // require vector and independent evaluator for kr (on faces!)
  Teuchos::ParameterList kre_list;
  kre_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  kre_list.setName("k_relative");
  auto kr_eval = Teuchos::rcp(new DiagIndependent(kre_list));
  S.SetEvaluator("k_relative", Tags::DEFAULT, kr_eval);

  // require boundary conditions
  Teuchos::ParameterList bce_list;
  bce_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  bce_list.setName("bcs");
  auto bc_eval = Teuchos::rcp(new BCsIndependent(bce_list));
  S.SetEvaluator("bcs", Tags::DEFAULT, bc_eval);

  // require the local operator and rhs
  Teuchos::ParameterList Ae_list;
  Ae_list.sublist("verbose object").set<std::string>("verbosity level", "high");
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
  S.SetEvaluator("A_local", Tags::DEFAULT, A_eval);
  S.SetEvaluator("A_rhs", Tags::DEFAULT, A_eval);

  // require secondary evaluator for r = b - Ax
  Teuchos::ParameterList re_list;
  re_list.sublist("verbose object").set<std::string>("verbosity level", "high");
  re_list.set("diagonal primary x key", "x");
  re_list.set("diagonal local operators keys", Teuchos::Array<std::string>(1, "A_local"));
  re_list.set("diagonal local operator rhss keys", Teuchos::Array<std::string>(1, "A_rhs"));
  re_list.set("additional rhss keys", Teuchos::Array<std::string>(1, "b"));
  re_list.set("rhs coefficients", Teuchos::Array<double>(1, -1.0));
  re_list.set("tag", "");

  re_list.set<std::string>("preconditioning method", "boomer amg");
  auto& amg_p = re_list.sublist("boomer amg parameters");
  amg_p.set("tolerance", 0.0);
  amg_p.set("verbosity", 0);

  re_list.set<std::string>("iterative method", "pcg");

  re_list.setName("residual");
  auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
  S.SetEvaluator("residual", Tags::DEFAULT, r_eval);
  S.Require<CompositeVector, CompositeVectorSpace>("residual", Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetGhosted(true);

  // require the derivative of r with respect to x.  The derivative of
  // OperatorApply WRT x is a linear operator.  We require assembly to invert,
  // and so the output of this is the assembled matrix.
  S.RequireDerivative<Operators::Operator, Operators::Operator_Factory>(
    "residual", Tags::DEFAULT, "x", Tags::DEFAULT);


  // run
  // -- setup and init
  S.Setup();
  S.GetW<CompositeVector>("x", Tags::DEFAULT, "x").PutScalar(1.);
  S.GetRecordW("x", Tags::DEFAULT, "x").set_initialized();
  S.Initialize();

  // -- update residual
  int updated = S.GetEvaluator("residual").Update(S, "pk");
  CHECK(updated);

  {
    // -- check error of the initial guess
    double error(0.);
    auto& r = S.Get<CompositeVector>("residual", Tags::DEFAULT);
    r.NormInf(&error);
    std::cout << "Error = " << error << std::endl;
  }

  // -- update derivative
  updated = S.GetEvaluator("residual").UpdateDerivative(S, "pk", "x", Tags::DEFAULT);
  CHECK(updated);

  {
    // -- get the derivative
    const auto& lin_op =
      S.GetDerivative<Operators::Operator>("residual", Tags::DEFAULT, "x", Tags::DEFAULT);

    // -- apply the inverse
    auto& r = S.Get<CompositeVector>("residual", Tags::DEFAULT);
    CompositeVector dx(r.Map());
    lin_op.ApplyInverse(r, dx);

    // -- update x
    S.GetW<CompositeVector>("x", Tags::DEFAULT, "x").Update(-1.0, dx, 1.0);
  }

  {
    // -- mark x as changed
    auto lx_eval = S.GetEvaluatorPtr("x", Tags::DEFAULT);
    auto x_primary_eval =
      Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(lx_eval);
    CHECK(x_primary_eval.get());
    x_primary_eval->SetChanged();
  }

  // -- update r again!
  updated = S.GetEvaluator("residual").Update(S, "pk");
  CHECK(updated);

  {
    // -- check error after the preconditioner is applied
    double error;
    auto& r = S.Get<CompositeVector>("residual", Tags::DEFAULT);
    r.NormInf(&error);
    std::cout << "Error = " << error << std::endl;
    CHECK_CLOSE(0.0, error, 1.e-3);
  }
}


SUITE(EVALUATOR_ON_OP)
{
  // Apply a non-diagonal operator, including boundary conditions
  TEST(OP_APPLY_DIFFUSION_FV)
  {
    test("fv: default");
  }
  // TEST(OP_APPLY_DIFFUSION_MFD) { test("mfd: two-point flux approximation"); }
  // TEST(OP_APPLY_DIFFUSION_NLFV) { test("nlfv: default"); }

  // Invert a non-diagonal operator, determining the solution to a linear problem
  TEST(OP_APPLY_DIFFUSION_FV_INVERSE)
  {
    test_inverse("fv: default");
  }
  // TEST(OP_APPLY_DIFFUSION_MFD_INVERSE) { test_inverse("mfd: two-point flux approximation"); }
  // TEST(OP_APPLY_DIFFUSION_NLFV_INVERSE) { test_inverse("nlfv: default"); }
}
