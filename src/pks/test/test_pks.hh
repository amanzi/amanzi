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

/*
ODE Test PKs

Test PKs:
  A: du/dt = 1     ( u = t )
  B: du/dt = u     ( u = exp(t) )
  C: du/dt = 2 * u * t  ( u = exp(t^2)  )

*/

#define DEBUG 1

using namespace Amanzi;

class DudtEvaluatorA
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  DudtEvaluatorA(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype(plist)
  {
    dependencies_.emplace_back(
      std::make_pair(Key("primaryA"), my_keys_[0].second));
  }

  DudtEvaluatorA(const DudtEvaluatorA& other) = default;
  Teuchos::RCP<Evaluator> Clone() const
  {
    return Teuchos::rcp(new DudtEvaluatorA(*this));
  }

 protected:
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
  {
    results[0]->putScalar(1.);
  }

  void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                  const Key& wrt_tag,
                                  const std::vector<CompositeVector*>& results)
  {
    AMANZI_ASSERT(std::make_pair(wrt_key, wrt_tag) == *dependencies_.begin());
    results[0]->putScalar(0.);
  }
};

class DudtEvaluatorB
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  DudtEvaluatorB(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype(plist)
  {
    dependencies_.emplace_back(
      std::make_pair(Key("primaryB"), my_keys_[0].second));
  }

  DudtEvaluatorB(const DudtEvaluatorB& other) = default;
  Teuchos::RCP<Evaluator> Clone() const
  {
    return Teuchos::rcp(new DudtEvaluatorB(*this));
  }

 protected:
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
  {
    results[0]->assign(S.Get<CompositeVector>(dependencies_.begin()->first,
            dependencies_.begin()->second));
  }

  void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                  const Key& wrt_tag,
                                  const std::vector<CompositeVector*>& results)
  {
    AMANZI_ASSERT(std::make_pair(wrt_key, wrt_tag) == *dependencies_.begin());
    results[0]->putScalar(1.);
  }
};

class DudtEvaluatorC
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  DudtEvaluatorC(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype(plist)
  {
    dependencies_.emplace_back(
      std::make_pair(Key("primaryC"), my_keys_[0].second));
    dependencies_.emplace_back(
      std::make_pair(Key("scaling"), my_keys_[0].second));
  }

  DudtEvaluatorC(const DudtEvaluatorC& other) = default;
  Teuchos::RCP<Evaluator> Clone() const
  {
    return Teuchos::rcp(new DudtEvaluatorC(*this));
  }

 protected:
  void Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
  {
    const auto& u = S.Get<CompositeVector>("primaryC", my_keys_[0].second);
    const auto& scaling = S.Get<CompositeVector>("scaling", my_keys_[0].second);
    results[0]->elementWiseMultiply(1.0, scaling, u, 0.);
  }

  void EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                                  const Key& wrt_tag,
                                  const std::vector<CompositeVector*>& results)
  {
    if (wrt_key == "primaryC") {
      results[0]->assign(S.Get<CompositeVector>("scaling", my_keys_[0].second));
    } else if (wrt_key == "scaling") {
      results[0]->assign(S.Get<CompositeVector>("primaryC", my_keys_[0].second));
    } else {
      AMANZI_ASSERT(false);
    }
  }
};

template <class Base_t, class Eval_t>
class PK_ODE_Explicit : public Base_t {
 public:
  PK_ODE_Explicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
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

  void Initialize()
  {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
      .putScalar(1.);
    this->S_->GetRecordW(this->key_, "", this->key_).set_initialized();
  }

  void FunctionalTimeDerivative(double t, const TreeVector& u, TreeVector& f)
  {
    this->S_->set_time(tag_inter_, t);
    // NOTE: Currently, we do not allow time integrators to introduce their own
    // tags -- time integrators have no knowledge of tags.  This makes it easy
    // for time integrators to break our state model, by simply
    // copy-constructing the input vector and giving us the copy instead of the
    // vector at the right tag.
    //    
    // For now, this ASSERT checks to ensure that the time integrator is giving
    // us what we expect.
    AMANZI_ASSERT(u.Data() ==
                  this->S_->template GetPtr<CompositeVector>(this->key_, this->tag_inter_));
    //
    // In the future, we should either have:
    //  1. The time integrator is aware of state, and works with tags.

    //  2. Call S_->SetPtr() and put the copied data into state.  This could be
    //     problematic, and would need some real testing, but might be the most
    //     convenient approach.
    //
    // For now, we also then need to set the value as changed, since something
    // may have happened at the inter tag.
    this->ChangedSolutionPK(tag_inter_);
    // END NOTE AND NONOBVIOUS CODE

    // evaluate the derivative
    this->S_->GetEvaluator(dudt_key_, tag_inter_).Update(*this->S_, this->name());
    *f.Data() = this->S_->template Get<CompositeVector>(dudt_key_, tag_inter_);

    std::cout << std::setprecision(16) << "  At time t = " << t
              << ": u = " << (u.Data()->ViewComponent<AmanziDefaultHost>("cell", false))(0, 0)
              << ", du/dt = " << (f.Data()->ViewComponent<AmanziDefaultHost>("cell", false))(0, 0)
              << std::endl;
  }

 protected:
  Key dudt_key_;
  using Base_t::tag_inter_;

  using Base_t::tag_new_;
  using Base_t::tag_old_;
};

template <class Base_t, class Eval_t>
class PK_ODE_Implicit : public Base_t {
 public:
  PK_ODE_Implicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
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

  void Initialize()
  {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
      .putScalar(1.);
    this->S_->GetRecordW(this->key_, this->key_).set_initialized();
  }

  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f)
  {
    AMANZI_ASSERT(std::abs(t_old - S_->time(tag_old_)) < 1.e-12);
    AMANZI_ASSERT(std::abs(t_new - S_->time(tag_new_)) < 1.e-12);

    // checking correct implementation of implicit time integrator (see note at
    // Implicit::FunctionalTimeDerivative)
    AMANZI_ASSERT(u_old->Data() ==
                  this->S_->template GetPtr<CompositeVector>(this->key_, tag_old_));
    AMANZI_ASSERT(u_new->Data() ==
                  this->S_->template GetPtr<CompositeVector>(this->key_, tag_new_));
    

    this->S_->GetEvaluator(dudt_key_, tag_new_).Update(*this->S_, this->name());
    f->Data()->assign(this->S_->template Get<CompositeVector>(dudt_key_, tag_new_));

    std::cout << "  At time t = " << t_new
              << ": u = " << (u_new->Data()->ViewComponent<AmanziDefaultHost>("cell", false))(0, 0)
              << ", du/dt = " << (f->Data()->ViewComponent<AmanziDefaultHost>("cell", false))(0, 0)
              << std::endl;

    double dt = t_new - t_old;
    f->update(1.0 / dt, *u_new, -1.0 / dt, *u_old, -1.0);
#if DEBUG
    {
      auto hv = f->Data()->ViewComponent<AmanziDefaultHost>("cell", false);
      std::cout << "R: res: " << hv(0,0)  << std::endl;
    }
#endif    
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
                          Teuchos::RCP<TreeVector> Pr)
  {
    Pr->Data()->assign(S_->template GetDerivative<CompositeVector>(
        this->dudt_key_, tag_new_, this->key_, tag_new_));
#if DEBUG
    {
      auto hv = r->Data()->ViewComponent<AmanziDefaultHost>("cell", false);
      std::cout << "PC: r: " << hv(0,0)  << std::endl;
    }

    {
      auto hv = Pr->Data()->ViewComponent<AmanziDefaultHost>("cell", false);
      std::cout << "PC: dRHS/du: " << hv(0,0)  << std::endl;
    }
#endif
  
    CompositeVector shift(r->Data()->getMap());
    shift.putScalar(1.0 / dt_);
    Pr->Data()->update(1., shift, 1.);

#if DEBUG
    {
      auto hv = Pr->Data()->ViewComponent<AmanziDefaultHost>("cell", false);
      std::cout << "PC: Jac: " << hv(0,0)  << std::endl;
    }
#endif

    //TreeVector tmp(Pr->getMap());
    Pr->reciprocal(*Pr);

#if DEBUG
    {
      auto hv = Pr->Data()->ViewComponent<AmanziDefaultHost>("cell", false);
      std::cout << "PC: 1/Jac: " << hv(0,0)  << std::endl;
    }
#endif

    Pr->elementWiseMultiply(1.0, *Pr, *r, 0.);

#if DEBUG
    {
      auto hv = Pr->Data()->ViewComponent<AmanziDefaultHost>("cell", false);
      std::cout << "PC: du: " << hv(0,0)  << std::endl;
    }
#endif
    return 0;
  }

  double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
  {
    AMANZI_ASSERT(u->Data() ==
                  this->S_->template GetPtr<CompositeVector>(this->key_, tag_new_));
    double norm = du->normInf();
    return norm;
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt)
  {
    AMANZI_ASSERT(u->Data() ==
                  this->S_->template GetPtr<CompositeVector>(this->key_, tag_new_));
    AMANZI_ASSERT(std::abs(this->S_->template Get<double>("time", tag_new_) -
                           this->S_->template Get<double>("time", tag_old_) - dt) < 1.e-10);
                  
    dt_ = dt;
    
    this->S_->GetEvaluator(dudt_key_, tag_new_)
      .UpdateDerivative(*this->S_, this->name(), this->key_, tag_new_);
  }

  void ChangedSolution() { this->ChangedSolutionPK(tag_new_); }

  void UpdateContinuationParameter(double lambda){};

 protected:
  double dt_;
  Key dudt_key_;

  using Base_t::S_;
  using Base_t::tag_new_;
  using Base_t::tag_old_;
};
