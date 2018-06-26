#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "EvaluatorAlgebraic.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Adaptors.hh"
#include "PK_Default.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"

/*
ODE Test PKs

Test PKs:
  A: du/dt = 1     ( u = t )
  B: du/dt = u     ( u = exp(t) )
  C: du/dt = 2 * u * t  ( u = exp(t^2)  )

*/

using namespace Amanzi;

class DudtEvaluatorA
    : public EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace> {
public:
  DudtEvaluatorA(Teuchos::ParameterList &plist) : EvaluatorAlgebraic(plist) {
    dependencies_.emplace_back(std::make_pair(Key("primaryA"), my_tag_));
  }

  DudtEvaluatorA(const DudtEvaluatorA &other) = default;
  Teuchos::RCP<Evaluator> Clone() const {
    return Teuchos::rcp(new DudtEvaluatorA(*this));
  }

protected:
  void Evaluate_(const State &S, CompositeVector &result) {
    result.PutScalar(1.);
  }

  void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                  const Key &wrt_tag, CompositeVector &result) {
    AMANZI_ASSERT(std::make_pair(wrt_key, wrt_tag) == *dependencies_.begin());
    result.PutScalar(0.);
  }
};

class DudtEvaluatorB
    : public EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace> {
public:
  DudtEvaluatorB(Teuchos::ParameterList &plist) : EvaluatorAlgebraic(plist) {
    dependencies_.emplace_back(std::make_pair(Key("primaryB"), my_tag_));
  }

  DudtEvaluatorB(const DudtEvaluatorB &other) = default;
  Teuchos::RCP<Evaluator> Clone() const {
    return Teuchos::rcp(new DudtEvaluatorB(*this));
  }

protected:
  void Evaluate_(const State &S, CompositeVector &result) {
    result = S.Get<CompositeVector>(dependencies_.begin()->first,
                                    dependencies_.begin()->second);
  }

  void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                  const Key &wrt_tag, CompositeVector &result) {
    AMANZI_ASSERT(std::make_pair(wrt_key, wrt_tag) == *dependencies_.begin());
    result.PutScalar(1.);
  }
};

class DudtEvaluatorC
    : public EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace> {
public:
  DudtEvaluatorC(Teuchos::ParameterList &plist) : EvaluatorAlgebraic(plist) {
    dependencies_.emplace_back(std::make_pair(Key("primaryC"), my_tag_));
    dependencies_.emplace_back(std::make_pair(Key("scaling"), my_tag_));
  }

  DudtEvaluatorC(const DudtEvaluatorC &other) = default;
  Teuchos::RCP<Evaluator> Clone() const {
    return Teuchos::rcp(new DudtEvaluatorC(*this));
  }

protected:
  void Evaluate_(const State &S, CompositeVector &result) {
    const auto &u = S.Get<CompositeVector>("primaryC", my_tag_);
    const auto &scaling = S.Get<CompositeVector>("scaling", my_tag_);
    result.Multiply(1.0, scaling, u, 0.);
  }

  void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                  const Key &wrt_tag, CompositeVector &result) {
    if (wrt_key == "primaryC") {
      result = S.Get<CompositeVector>("scaling", my_tag_);
    } else if (wrt_key == "scaling") {
      result = S.Get<CompositeVector>("primaryC", my_tag_);
    } else {
      AMANZI_ASSERT(false);
    }
  }
};

template <class Base_t, class Eval_t> class PK_ODE_Explicit : public Base_t {
public:
  PK_ODE_Explicit(const Teuchos::RCP<Teuchos::ParameterList> &pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList> &global_plist,
                  const Teuchos::RCP<State> &S)
      : Base_t(pk_tree, global_plist, S) {
    dudt_key_ = this->key_ + "_t";
  }

  void Setup() {
    Base_t::Setup();

    this->S_
        ->template Require<CompositeVector, CompositeVectorSpace>(
            this->key_, tag_old_, this->key_)
        .SetMesh(this->mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    this->S_
        ->template Require<CompositeVector, CompositeVectorSpace>(
            dudt_key_, this->tag_inter_, dudt_key_)
        .SetMesh(this->mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList plist(dudt_key_);
    plist.set("tag", this->tag_inter_);
    auto dudt_eval = Teuchos::rcp(new Eval_t(plist));
    this->S_->SetEvaluator(dudt_key_, this->tag_inter_, dudt_eval);
  }

  void Initialize() {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
        .PutScalar(1.);
    this->S_->GetRecordW(this->key_, "", this->key_).set_initialized();
  }

  void Functional(double t, const TreeVector &u, TreeVector &f) {
    // these calls are here because the explicit ti is not aware that it should
    // call it
    this->S_->set_time(tag_inter_, t);
    this->ChangedSolutionPK(tag_inter_);

    // evaluate the derivative
    this->S_->GetEvaluator(dudt_key_, tag_inter_)
        .Update(*this->S_, this->name());
    *f.Data() = this->S_->template Get<CompositeVector>(dudt_key_, tag_inter_);

    std::cout << std::setprecision(16) << "  At time t = " << t
              << ": u = " << (*u.Data()->ViewComponent("cell", false))[0][0]
              << ", du/dt = " << (*f.Data()->ViewComponent("cell", false))[0][0]
              << std::endl;
  }

protected:
  Key dudt_key_;
  using Base_t::tag_inter_;

  using Base_t::tag_new_;
  using Base_t::tag_old_;
};

template <class Base_t, class Eval_t> class PK_ODE_Implicit : public Base_t {
public:
  PK_ODE_Implicit(const Teuchos::RCP<Teuchos::ParameterList> &pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList> &global_plist,
                  const Teuchos::RCP<State> &S)
      : Base_t(pk_tree, global_plist, S) {
    dudt_key_ = this->key_ + "_t";
  }

  void Setup() {
    Base_t::Setup();

    this->S_
        ->template Require<CompositeVector, CompositeVectorSpace>(
            this->key_, tag_new_, this->key_)
        .SetMesh(this->mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    this->S_
        ->template Require<CompositeVector, CompositeVectorSpace>(
            dudt_key_, tag_new_, dudt_key_)
        .SetMesh(this->mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    this->S_->template RequireDerivative<CompositeVector, CompositeVectorSpace>(
        dudt_key_, tag_new_, this->key_, tag_new_);

    Teuchos::ParameterList plist(dudt_key_);
    plist.set("tag", tag_new_);
    auto dudt_eval = Teuchos::rcp(new Eval_t(plist));
    this->S_->SetEvaluator(dudt_key_, tag_new_, dudt_eval);
  }

  void Initialize() {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
        .PutScalar(1.);
    this->S_->GetRecordW(this->key_, this->key_).set_initialized();
  }

  void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                  Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
    AMANZI_ASSERT(std::abs(t_old - S_->time(tag_old_)) < 1.e-12);
    AMANZI_ASSERT(std::abs(t_new - S_->time(tag_new_)) < 1.e-12);

    this->S_->GetEvaluator(dudt_key_, tag_new_)
        .Update(*this->S_, this->name());
    *f->Data() = this->S_->template Get<CompositeVector>(dudt_key_, tag_new_);

    std::cout << "  At time t = " << t_new << ": u = "
              << (*u_new->Data()->ViewComponent("cell", false))[0][0]
              << ", du/dt = "
              << (*f->Data()->ViewComponent("cell", false))[0][0] << std::endl;

    double dt = t_new - t_old;
    f->Update(1.0 / dt, *u_new, -1.0 / dt, *u_old, -1.0);
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                          Teuchos::RCP<TreeVector> Pu) {
    *Pu->Data() = S_->template GetDerivative<CompositeVector>(
        this->dudt_key_, tag_new_, this->key_, tag_new_);
    Pu->Data()->Shift(1.0 / h_);
    Pu->ReciprocalMultiply(1.0, *Pu, *u, 0.);
    return 0;
  }

  double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<const TreeVector> du) {
    double norm = 0.;
    du->NormInf(&norm);
    std::cout << "     error = " << norm << std::endl;
    return norm;
  }

  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up,
                            double h) {
    h_ = h;
    this->S_->GetEvaluator(dudt_key_, tag_new_)
        .UpdateDerivative(*this->S_, this->name(), this->key_, tag_new_);
  }

  void ChangedSolution() { this->ChangedSolutionPK(tag_new_); }

  void UpdateContinuationParameter(double lambda){};

protected:
  double h_;
  Key dudt_key_;

  using Base_t::S_;
  using Base_t::tag_new_;
  using Base_t::tag_old_;
};
