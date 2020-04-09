/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

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
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorIndependent.hh"
#include "evaluator/EvaluatorPrimary.hh"
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

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new BIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    double cv = s.GetMesh()->cell_volume(0, false);
    auto& b = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    b.ViewComponent("cell", false)->putScalar(-4. * cv);
  }
};

class XIndependent
  : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  using EvaluatorIndependent<CompositeVector,
                             CompositeVectorSpace>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new XIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    auto& x = s.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    auto& x_c = *x.ViewComponent("cell", false);
    for (int c = 0; c != x_c.getLocalLength(); ++c) {
      AmanziGeometry::Point cc = x.Mesh()->cell_centroid(c);
      x_c[0][c] = cc[0] * cc[0] + cc[1] * cc[1]; // x^2 + y^2
    }

    if (x.HasComponent("face")) {
      auto& x_f = *x.ViewComponent("face", false);
      for (int f = 0; f != x_f.getLocalLength(); ++f) {
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

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new DiagIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    s.GetW<CompositeVector>(my_key_, my_tag_, my_key_).putScalar(1.);
  }
};

class KIndependent
  : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  using EvaluatorIndependent<TensorVector,
                             TensorVector_Factory>::EvaluatorIndependent;

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
      k.putScalar(1.0);
    }
  }
};

class BCsIndependent
  : public EvaluatorIndependent<Operators::BCs, Operators::BCs_Factory> {
 public:
  using EvaluatorIndependent<Operators::BCs,
                             Operators::BCs_Factory>::EvaluatorIndependent;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new BCsIndependent(*this));
  };

 protected:
  virtual void Update_(State& s) override
  {
    auto& bcs = s.GetW<Operators::BCs>(my_key_, my_tag_, my_key_);

    auto mesh = bcs.mesh();
    auto kind = bcs.kind();

    auto& model = bcs.bc_model();
    auto& value = bcs.bc_value();

    for (auto& bc : model) bc = Operators::OPERATOR_BC_NONE;
    for (auto& val : value) val = 0.;

    // set all exterior faces to dirichlet 0
    int nfaces_owned =
      mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    AmanziMesh::Entity_ID_List cells;
    for (int f = 0; f != nfaces_owned; ++f) {
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      if (cells.size() == 1) {
        model[f] = Operators::OPERATOR_BC_DIRICHLET;
        auto fc = mesh->face_centroid(f);
        value[f] = fc[0] * fc[0] + fc[1] * fc[1]; // x^2 + y^2
      }
    }
  }
};


void
test(const std::string& discretization)
{
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
  Teuchos::ParameterList re_list("residual");
  re_list.sublist("verbose object")
    .set<std::string>("verbosity level", "extreme");
  re_list.set("diagonal primary x key", "x");
  re_list.set("diagonal local operators keys",
              Teuchos::Array<std::string>(1, "A_local"));
  re_list.set("diagonal local operator rhss keys",
              Teuchos::Array<std::string>(1, "A_rhs"));
  re_list.set("additional rhss keys", Teuchos::Array<std::string>(1, "b"));
  re_list.set("rhs coefficients", Teuchos::Array<double>(1, -1.0));
  re_list.set("tag", "");
  auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
  S.SetEvaluator("residual", r_eval);
  S.Require<CompositeVector, CompositeVectorSpace>("residual", "")
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
  double error(0.);
  auto& r = S.Get<CompositeVector>("residual", "");
  error = r.normInf();
  std::cout << "Error = " << error << std::endl;
  CHECK(error != 0.0);
  CHECK_CLOSE(0.0, error, 1.e-3);
}


void
test_inverse(const std::string& discretization)
{
  auto comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  MeshFactory meshfac(comm);
  //  auto mesh = meshfac(-1.0, -1.0, 1.0, 1.0, 128, 128);
  auto mesh = meshfac(-1.0, -1.0, 1.0, 1.0, 4, 4);

  State S;
  S.RegisterDomainMesh(mesh);

  Teuchos::ParameterList es_list;
  es_list.sublist("verbose object")
    .set<std::string>("verbosity level", "extreme");
  es_list.setName("my_op");

  // NOTE: still need to require evaluators because factory is not tested
  // yet.  Just testing remove of DATA needs.

  // require primary evaluator for x
  Teuchos::ParameterList xe_list;
  xe_list.sublist("verbose object")
    .set<std::string>("verbosity level", "extreme");
  xe_list.setName("x");
  auto x_eval = Teuchos::rcp(
    new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(xe_list));
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

  // require secondary evaluator for r = b - Ax
  Teuchos::ParameterList re_list;
  re_list.sublist("verbose object")
    .set<std::string>("verbosity level", "extreme");
  re_list.set("diagonal primary x key", "x");
  re_list.set("diagonal local operators keys",
              Teuchos::Array<std::string>(1, "A_local"));
  re_list.set("diagonal local operator rhss keys",
              Teuchos::Array<std::string>(1, "A_rhs"));
  re_list.set("additional rhss keys", Teuchos::Array<std::string>(1, "b"));
  re_list.set("rhs coefficients", Teuchos::Array<double>(1, -1.0));
  re_list.sublist("preconditioner")
    .set<std::string>("preconditioner type", "boomer amg");
  re_list.set("tag", "");
  auto& amg_p =
    re_list.sublist("preconditioner").sublist("boomer amg parameters");
  amg_p.set("tolerance", 0.0);
  amg_p.set("verbosity", 3);

  re_list.setName("residual");
  auto r_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
  S.SetEvaluator("residual", r_eval);
  S.Require<CompositeVector, CompositeVectorSpace>("residual", "")
    .SetMesh(mesh)
    ->SetGhosted(true);

  // require the derivative of r with respect to x.  The derivative of
  // OperatorApply WRT x is a linear operator.  We require assembly to invert,
  // and so the output of this is the assembled matrix.
  S.RequireDerivative<Operators::Operator, Operators::Operator_Factory>(
    "residual", "", "x", "");


  // run
  // -- setup and init
  S.Setup();
  S.GetW<CompositeVector>("x", "", "x").putScalar(1.);
  S.GetRecordW("x", "", "x").set_initialized();
  S.Initialize();

  // -- update residual
  int updated = S.GetEvaluator("residual").Update(S, "pk");
  CHECK(updated);

  {
    // -- check error of the initial guess
    double error(0.);
    auto& r = S.Get<CompositeVector>("residual", "");
    error = r.normInf();
    std::cout << "Error = " << error << std::endl;
  }

  // -- update derivative
  updated = S.GetEvaluator("residual").UpdateDerivative(S, "pk", "x", "");
  CHECK(updated);

  {
    // -- get the derivative
    const auto& lin_op =
      S.GetDerivative<Operators::Operator>("residual", "", "x", "");

    // -- apply the inverse
    auto& r = S.Get<CompositeVector>("residual", "");
    CompositeVector dx(r.getMap());
    std::cout << "R:" << std::endl;
    r.Print(std::cout);
    lin_op.applyInverse(r, dx);
    std::cout << "dx:" << std::endl;
    dx.Print(std::cout);

    // -- update x
    S.GetW<CompositeVector>("x", "", "x").Update(-1.0, dx, 1.0);
    std::cout << "x corrected:" << std::endl;
    S.Get<CompositeVector>("x", "").Print(std::cout);
  }

  {
    // -- mark x as changed
    auto lx_eval = S.GetEvaluatorPtr("x", "");
    auto x_primary_eval = Teuchos::rcp_dynamic_cast<
      EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(lx_eval);
    CHECK(x_primary_eval.get());
    x_primary_eval->SetChanged();
  }

  // -- update r again!
  updated = S.GetEvaluator("residual").Update(S, "pk");
  CHECK(updated);

  {
    // -- check error after the preconditioner is applied
    double error(0.);
    auto& r = S.Get<CompositeVector>("residual", "");
    r.Print(std::cout);
    error = r.normInf();
    std::cout << "Error = " << error << std::endl;
    CHECK(error != 0.0);
    CHECK_CLOSE(0.0, error, 1.e-3);
  }
}


SUITE(EVALUATOR_ON_OP)
{
  // Apply a non-diagonal operator, including boundary conditions
  TEST(OP_APPLY_DIFFUSION_FV) { test("fv: default"); }
  TEST(OP_APPLY_DIFFUSION_MFD) { test("mfd: two-point flux approximation"); }
  //  TEST(OP_APPLY_DIFFUSION_NLFV) { test("nlfv: default"); }

  // Invert a non-diagonal operator, determining the solution to a linear
  // problem
  TEST(OP_APPLY_DIFFUSION_FV_INVERSE) { test_inverse("fv: default"); }
  TEST(OP_APPLY_DIFFUSION_MFD_INVERSE)
  {
    test_inverse("mfd: two-point flux approximation");
  }
  //  TEST(OP_APPLY_DIFFUSION_NLFV_INVERSE) { test_inverse("nlfv: default"); }
}
