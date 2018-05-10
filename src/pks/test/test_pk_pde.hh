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

#include "Evaluator_OperatorApply.hh"
#include "Operator_Factory.hh"

/*
PDE Test PKs

1. Diffusion equation: du/dt - div k grad u = q

*/

using namespace Amanzi;


template <class Base_t> class PK_PDE_Explicit : public Base_t {
public:
  PK_PDE_Explicit(const Teuchos::RCP<Teuchos::ParameterList> &pk_tree,
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
        .SetMesh(this->mesh_);

    this->S_
        ->template Require<CompositeVector, CompositeVectorSpace>(
            dudt_key_, this->tag_inter_, dudt_key_)
        .SetMesh(this->mesh_);

    Teuchos::ParameterList re_list(dudt_key_);
    re_list.sublist("verbose object")
        .set<std::string>("verbosity level", "high");
    re_list.set("diagonal primary x key", this->key_);
    re_list.set("diagonal local operators keys", Teuchos::Array<std::string>(1,"A_local"));
    re_list.set("diagonal local operator rhss keys", Teuchos::Array<std::string>(1,"A_rhs"));
    re_list.set("additional rhss keys", Teuchos::Array<std::string>(1,"source_cv"));
    re_list.set("rhs coefficients", Teuchos::Array<double>(1,1.0));
    re_list.set("tag", this->tag_inter_);
    auto dudt_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
    this->S_->SetEvaluator(dudt_key_, this->tag_inter_, dudt_eval);
  }

  void Initialize() {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
        .PutScalar(0.);
    this->S_->GetRecordW(this->key_, "", this->key_).set_initialized();
  }

  void Dudt(double t, const TreeVector &u, TreeVector &f) {
    // these calls are here because the explicit ti is not aware that it should
    // call it
    this->S_->set_time(tag_inter_, t);
    this->ChangedSolutionPK(tag_inter_);

    // evaluate the derivative
    this->S_->GetEvaluator(dudt_key_, tag_inter_)
        .Update(*this->S_, this->name());
    this->S_->GetEvaluator("cell_volume", tag_inter_)
        .Update(*this->S_, this->name());
    f.Data()->ReciprocalMultiply(1.0, this->S_->template Get<CompositeVector>("cell_volume", tag_inter_),
            this->S_->template Get<CompositeVector>(dudt_key_, tag_inter_), 0.);            

    std::cout << std::setprecision(16) << "  At time t = " << t
              << ": u = " << (*u.Data()->ViewComponent("cell", false))[0][60]
              << ", du/dt = " << (*f.Data()->ViewComponent("cell", false))[0][60]
              << ", q = " << (*this->S_->template Get<CompositeVector>("source",tag_inter_).ViewComponent("cell",false))[0][60]
              << std::endl;
  }

protected:
  Key dudt_key_;
  using Base_t::tag_inter_;

  using Base_t::tag_new_;
  using Base_t::tag_old_;
};



template <class Base_t> class PK_PDE_Implicit : public Base_t {
public:
  PK_PDE_Implicit(const Teuchos::RCP<Teuchos::ParameterList> &pk_tree,
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
        .SetMesh(this->mesh_);

    this->S_
        ->template Require<CompositeVector, CompositeVectorSpace>(
            dudt_key_, tag_new_, dudt_key_)
        .SetMesh(this->mesh_);

    this->S_->template RequireDerivative<Operators::Operator,Operators::Operator_Factory>(
        dudt_key_, tag_new_, this->key_, tag_new_);


    Teuchos::ParameterList re_list(dudt_key_);
    re_list.sublist("verbose object")
        .set<std::string>("verbosity level", "high");
    re_list.set("diagonal primary x key", this->key_);
    re_list.set("diagonal local operators keys", Teuchos::Array<std::string>(1,"A_local"));
    re_list.set("diagonal local operator rhss keys", Teuchos::Array<std::string>(1,"A_rhs"));
    re_list.set("additional rhss keys", Teuchos::Array<std::string>(1,"source_cv"));
    re_list.set("rhs coefficients", Teuchos::Array<double>(1,1.0));
    re_list.set("tag", this->tag_new_);
    auto dudt_eval = Teuchos::rcp(new Evaluator_OperatorApply(re_list));
    this->S_->SetEvaluator(dudt_key_, this->tag_new_, dudt_eval);
    //    this->S_->template RequireEvaluator(dudt_key_, tag_new_);
  }

  void Initialize() {
    Base_t::Initialize();
    this->S_->template GetW<CompositeVector>(this->key_, "", this->key_)
        .PutScalar(0.);
    this->S_->GetRecordW(this->key_, "", this->key_).set_initialized();
  }



  void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                  Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
    ASSERT(std::abs(t_old - this->S_->time(tag_old_)) < 1.e-12);
    ASSERT(std::abs(t_new - this->S_->time(tag_new_)) < 1.e-12);

    this->S_->GetEvaluator(dudt_key_, tag_new_)
        .Update(*this->S_, this->name());
    *f->Data() = this->S_->template Get<CompositeVector>(dudt_key_, tag_new_);

    std::cout << "  At time t = " << t_new << ": u = "
              << (*u_new->Data()->ViewComponent("cell", false))[0][60]
              << ", du/dt = "
              << (*f->Data()->ViewComponent("cell", false))[0][60] << std::endl;

    double dt = t_new - t_old;
    f->Update(1.0 / dt, *u_new, -1.0 / dt, *u_old, -1.0);
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                          Teuchos::RCP<TreeVector> Pu) {
    const auto& lin_op = this->S_->template GetDerivative<Operators::Operator>(dudt_key_, tag_new_, this->key_, tag_new_);
    lin_op.ApplyInverse(*u->Data(), *Pu->Data());
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
    this->S_->GetEvaluator(dudt_key_, tag_new_)
        .UpdateDerivative(*this->S_, this->name(), this->key_, tag_new_);
  }

  void ChangedSolution() { this->ChangedSolutionPK(tag_new_); }

  void UpdateContinuationParameter(double lambda){};
  
  // void Dudt(double t, const TreeVector &u, TreeVector &f) {
  //   // these calls are here because the explicit ti is not aware that it should
  //   // call it
  //   this->S_->set_time(tag_inter_, t);
  //   this->ChangedSolutionPK(tag_inter_);

  //   // evaluate the derivative
  //   this->S_->GetEvaluator(dudt_key_, tag_inter_)
  //       .Update(*this->S_, this->name());
  //   this->S_->GetEvaluator("cell_volume", tag_inter_)
  //       .Update(*this->S_, this->name());
  //   f.Data()->ReciprocalMultiply(1.0, this->S_->template Get<CompositeVector>("cell_volume", tag_inter_),
  //           this->S_->template Get<CompositeVector>(dudt_key_, tag_inter_), 0.);            

  //   std::cout << std::setprecision(16) << "  At time t = " << t
  //             << ": u = " << (*u.Data()->ViewComponent("cell", false))[0][60]
  //             << ", du/dt = " << (*f.Data()->ViewComponent("cell", false))[0][60]
  //             << ", q = " << (*this->S_->template Get<CompositeVector>("source",tag_inter_).ViewComponent("cell",false))[0][60]
  //             << std::endl;
  // }

protected:
  Key dudt_key_;

  using Base_t::tag_new_;
  using Base_t::tag_old_;
};
