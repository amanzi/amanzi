/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Adaptors.hh"
#include "PK_Default.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"

#include "Evaluator_OperatorApply.hh"
#include "Evaluator_PDE_Accumulation.hh"
#include "Operator_Factory.hh"


static const double PI_2 = 1.5707963267948966;

/*
PDE Test PKs

1. Diffusion equation: du/dt - div k grad u = q

*/

using namespace Amanzi;


template <class Base_t>
class PK_PDE_Explicit : public Base_t {
 public:
  PK_PDE_Explicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                  const Teuchos::RCP<State>& S)
    : Base_t(pk_tree, global_plist, S)
  {
    dudt_key_ = this->key_ + "_t";
  }

  void Setup()
  {
    Base_t::Setup();

    this->S_
      ->template Require<CompositeVector, CompositeVectorSpace>(
        this->key_, tag_old_, this->key_)
      .SetMesh(this->mesh_);

    this->S_
      ->template Require<CompositeVector, CompositeVectorSpace>(
        dudt_key_, this->tag_inter_, dudt_key_)
      .SetMesh(this->mesh_);

    Teuchos::ParameterList re_list(dudt_key_);
    re_list.sublist("verbose object")
      .set<std::string>("verbosity level", "high");
    re_list.set("tag", this->tag_inter_);
    re_list.set("diagonal primary x key", this->key_);
    re_list.set("diagonal local operators keys",
                Teuchos::Array<std::string>(1, "A_local"));
    re_list.set("diagonal local operator rhss keys",
                Teuchos::Array<std::string>(1, "A_rhs"));
    re_list.set("additional rhss keys",
                Teuchos::Array<std::string>(1, "source_cv"));
    re_list.set("rhs coefficients", Teuchos::Array<double>(1, -1.0));
    auto dudt_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
    this->S_->SetEvaluator(this->dudt_key_, this->tag_inter_, dudt_eval);
  }

  void Initialize()
  {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
      .putScalar(0.);
    this->S_->GetRecordW(this->key_, "", this->key_).set_initialized();
  }

  void FunctionalTimeDerivative(double t, const TreeVector& u, TreeVector& f)
  {
    // these calls are here because the explicit ti is not aware that it should
    // call it
    this->S_->set_time(tag_inter_, t);
    this->ChangedSolutionPK(tag_inter_);

    // evaluate the derivative
    this->S_->GetEvaluator(dudt_key_, tag_inter_)
      .Update(*this->S_, this->name());
    const auto& dudt =
      this->S_->template Get<CompositeVector>(dudt_key_, tag_inter_);

    // evaluate cell volume
    this->S_->GetEvaluator("cell_volume", tag_inter_)
      .Update(*this->S_, this->name());
    const auto& cv =
      this->S_->template Get<CompositeVector>("cell_volume", tag_inter_);

    // dudt = -dudt_eval / cv
    f.Data()->reciprocal(cv);
    f.Data()->elementWiseMultiply(-1.0, dudt, *f.Data(), 0.);

    std::cout << "Dudt at t = " << t << std::endl;
    std::cout << "================================================"
              << std::endl;
    std::cout << "  u = " << std::endl;
    u.Print(std::cout);
  }

 protected:
  Key dudt_key_;
  using Base_t::tag_inter_;

  using Base_t::tag_new_;
  using Base_t::tag_old_;
};


template <class Base_t>
class PK_PDE_Implicit : public Base_t {
 public:
  PK_PDE_Implicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                  const Teuchos::RCP<State>& S)
    : Base_t(pk_tree, global_plist, S)
  {
    res_key_ = this->key_ + "_res";
    conserved_key_ =
      this->plist_->template get<std::string>("conserved quantity", "u");
  }

  void Setup()
  {
    Base_t::Setup();

    this->S_
      ->template Require<CompositeVector, CompositeVectorSpace>(
        this->key_, tag_new_, this->key_)
      .SetMesh(this->mesh_);

    this->S_
      ->template Require<CompositeVector, CompositeVectorSpace>(
        res_key_, tag_new_, res_key_)
      .SetMesh(this->mesh_);

    this->S_->template RequireDerivative<Operators::Operator,
                                         Operators::Operator_Factory>(
      res_key_, tag_new_, this->key_, tag_new_);

    // set up the accumulation
    Teuchos::ParameterList acc_list("accumulation");
    acc_list.sublist("verbose object")
      .set<std::string>("verbosity level", "high");
    acc_list.set("conserved quantity key", this->key_);
    acc_list.set("tag", this->tag_new_);
    acc_list.set("tag old", this->tag_old_);
    acc_list.set("tag new", this->tag_new_);
    auto acc_eval = Teuchos::rcp(new Evaluator_PDE_Accumulation(acc_list));
    this->S_->SetEvaluator("accumulation", this->tag_new_, acc_eval);

    // set up the residual evaluator
    Teuchos::ParameterList re_list(res_key_);
    re_list.sublist("verbose object")
      .set<std::string>("verbosity level", "high");
    re_list.set("tag", this->tag_new_);

    // -- diffusion operator
    re_list.set("diagonal primary x key", this->key_);
    re_list.set("diagonal local operators keys",
                Teuchos::Array<std::string>(1, "A_local"));
    re_list.set("diagonal local operator rhss keys",
                Teuchos::Array<std::string>(1, "A_rhs"));

    // -- rhs for source and accumulation term
    auto rhss = Teuchos::Array<std::string>(
      std::vector<std::string>{ "source_cv", "accumulation" });
    re_list.set("additional rhss keys", rhss);
    auto rhs_coefs = Teuchos::Array<double>(std::vector<double>{ -1.0, 1.0 });
    re_list.set("rhs coefficients", rhs_coefs);

    // -- preconditioner for evaluating inverse
    re_list.sublist("preconditioner").set("preconditioner type", "diagonal");
    /*
    re_list.sublist("preconditioner").set("preconditioner type", "boomer amg");
    re_list.sublist("preconditioner")
      .sublist("boomer amg parameters")
      .set("cycle iterations", 2);
    re_list.sublist("preconditioner")
      .sublist("boomer amg parameters")
      .set("smoother sweeps", 5);
    re_list.sublist("preconditioner")
      .sublist("boomer amg parameters")
      .set("tolerance", 0.0);
    re_list.sublist("preconditioner")
      .sublist("boomer amg parameters")
      .set("verbosity", 0);
    */
    auto res_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
    this->S_->SetEvaluator(res_key_, this->tag_new_, res_eval);
  }

  void Initialize()
  {
    Base_t::Initialize();
    auto& ic_cv =
      this->S_->template GetW<CompositeVector>(this->key_, "", this->key_);
    ic_cv.putScalar(0.);
    this->S_->GetRecordW(this->key_, "", this->key_).set_initialized();
  }

  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f)
  {
    AMANZI_ASSERT(std::abs(t_old - this->S_->time(tag_old_)) < 1.e-12);
    AMANZI_ASSERT(std::abs(t_new - this->S_->time(tag_new_)) < 1.e-12);

    std::cout << "Updating residual at time t = " << t_new << std::endl;
    std::cout << "================================================"
              << std::endl;
    // std::cout << "  u_old = ";
    // u_old->Print(std::cout);
    // std::cout << "  u_new = ";
    // u_new->Print(std::cout);
    this->S_->GetEvaluator(res_key_, tag_new_).Update(*this->S_, this->name());
    *f->Data() = this->S_->template Get<CompositeVector>(res_key_, tag_new_);

    // std::cout << "  res = ";
    // f->Print(std::cout);
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                          Teuchos::RCP<TreeVector> Pu)
  {
    const auto& lin_op = this->S_->template GetDerivative<Operators::Operator>(
      res_key_, tag_new_, this->key_, tag_new_);
    lin_op.applyInverse(*u->Data(), *Pu->Data());
    return 0;
  }

  double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
  {
    double norm = du->normInf();
    std::cout << "     error = " << norm << std::endl;
    return norm;
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
  {
    this->S_->GetEvaluator(res_key_, tag_new_)
      .UpdateDerivative(*this->S_, this->name(), this->key_, tag_new_);
  }

  void ChangedSolution() { this->ChangedSolutionPK(tag_new_); }

  void UpdateContinuationParameter(double lambda){};

 protected:
  Key res_key_;
  Key conserved_key_;

  using Base_t::tag_new_;
  using Base_t::tag_old_;
};
