#ifndef AMANZI_NOX_SOLVER_FN_BASE_
#define AMANZI_NOX_SOLVER_FN_BASE_

#include "NOX_Abstract_Group.H"
#include "NoxVector.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class VectorClass>
class AmanziGroup : public NOX::Abstract::Group {

 public:
  AmanziGroup(const Teuchos::RCP<SolverFnBase<VectorClass> > fn) {
    fn_ = fn;
    is_f_ = false;
    is_jac_ = false;
    is_dx_ = false;
  }

  NOX::Abstract::Group&
  operator=(const NOX::Abstract::Group& source) override {
    const AmanziGroup<VectorClass>& sourceGroup = dynamic_cast<const AmanziGroup<VectorClass>&>(source);

    fn_ = sourceGroup.fn_;
    if (Teuchos::nonnull(sourceGroup.x_))
      x_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(sourceGroup.x_->clone());
    if (Teuchos::nonnull(sourceGroup.f_))
      f_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(sourceGroup.f_->clone());
    if (Teuchos::nonnull(sourceGroup.dx_))
      dx_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(sourceGroup.dx_->clone());
    is_f_ = sourceGroup.is_f_;
    is_jac_ = sourceGroup.is_jac_;
    is_dx_ = sourceGroup.is_dx_;

    return *this;
  }

  void setX(const NOX::Abstract::Vector &y) override {
    fn_->ChangedSolution();
    x_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(y.clone());
    is_f_ = false;
    is_jac_ = false;
    is_dx_ = false;
  }

  // Compute x = grp.x + step * d
  void computeX(const NOX::Abstract::Group &grp, const NOX::Abstract::Vector &d, double step) override
  {
    const AmanziGroup<VectorClass>& sourceGroup = dynamic_cast<const AmanziGroup<VectorClass>&>(grp);
    assert(Teuchos::nonnull(x_));
    assert(Teuchos::nonnull(sourceGroup.x_));

    x_->update(1.,*sourceGroup.x_,step,d);

    is_f_ = false;
    is_jac_ = false;
    is_dx_ = false;
  }

  NOX::Abstract::Group::ReturnType
  computeF() override {
    assert(Teuchos::nonnull(x_));
    if (Teuchos::is_null(f_)) {
      f_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(x_->clone(NOX::ShapeCopy));
    }
    assert(Teuchos::nonnull(f_));
    fn_->Residual(x_->get_vector(), f_->get_vector());
    is_f_ = true;
    return NOX::Abstract::Group::Ok;
  }

  NOX::Abstract::Group::ReturnType
  computeJacobian() override {
    assert(Teuchos::nonnull(x_));
    fn_->UpdatePreconditioner(x_->get_vector());
    is_jac_ = true;
    return NOX::Abstract::Group::Ok;
  }

  // NOTE: Room for improvement for things like Eisenstat-Walker for tuning
  // tols by passing these params into a PK's linear operator
  NOX::Abstract::Group::ReturnType
  computeNewton(Teuchos::ParameterList& params) override {
    assert(Teuchos::nonnull(f_));
    if (Teuchos::is_null(dx_)) {
      dx_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(x_->clone(NOX::ShapeCopy));
    }
    assert(Teuchos::nonnull(dx_));
    assert(isJacobian());
    
    int ierr = fn_->ApplyPreconditioner(f_->get_vector(), dx_->get_vector());
    dx_->scale(-1.);
    if (ierr > 0) {
      return NOX::Abstract::Group::Ok;
    } else {
      return NOX::Abstract::Group::NotConverged;
    }
  }

  
  NOX::Abstract::Group::ReturnType
  applyJacobianInverse(Teuchos::ParameterList& params,
                       const NOX::Abstract::Vector& input,
                       NOX::Abstract::Vector& result) const override {
    assert(0);
    const NoxVector<VectorClass>* input_a =
        dynamic_cast<const NoxVector<VectorClass>*>(&input);
    NoxVector<VectorClass>* result_a =
        dynamic_cast<NoxVector<VectorClass>*>(&result);

    int ierr = fn_->ApplyPreconditioner(input_a->get_vector(), result_a->get_vector());
    dx_->scale(-1.);
    if (ierr > 0) {
      return NOX::Abstract::Group::Ok;
    } else {
      return NOX::Abstract::Group::NotConverged;
    }
  }
  
  bool isF() const override { return is_f_; }
  bool isJacobian() const override { return is_jac_; }

  const NOX::Abstract::Vector& getX() const override { assert(Teuchos::nonnull(x_)); return *x_; }
  const NOX::Abstract::Vector& getF() const override { assert(Teuchos::nonnull(f_)); return *f_; }
  double getNormF() const override { assert(Teuchos::nonnull(f_)); return f_->norm(); }
  const NOX::Abstract::Vector& getGradient() const override { assert(0); return *dx_; }
  const NOX::Abstract::Vector& getNewton() const override { assert(Teuchos::nonnull(dx_)); return *dx_; }
  Teuchos::RCP<NOX::Abstract::Group> clone(NOX::CopyType type=NOX::DeepCopy) const override
  {
    Teuchos::RCP<AmanziGroup<VectorClass> > cloneGroup = Teuchos::rcp(new AmanziGroup<VectorClass>(fn_));
    cloneGroup->fn_ = fn_;
    if (Teuchos::nonnull(x_))
      cloneGroup->x_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(x_->clone());
    if (Teuchos::nonnull(f_))
      cloneGroup->f_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(f_->clone());
    if (Teuchos::nonnull(dx_))
      cloneGroup->dx_ = Teuchos::rcp_dynamic_cast<NoxVector<VectorClass> >(dx_->clone());
    cloneGroup->is_f_ = is_f_;
    cloneGroup->is_jac_ = is_jac_;
    cloneGroup->is_dx_ = is_dx_;

    return cloneGroup;
  }

  // remaining defualt interfaces
  virtual Teuchos::RCP<const NOX::Abstract::Vector> getXPtr() const override { return x_; }
  virtual Teuchos::RCP<const NOX::Abstract::Vector> getFPtr() const override { return f_; }
  virtual Teuchos::RCP<const NOX::Abstract::Vector> getGradientPtr() const override { assert(0); return Teuchos::null; }
  virtual Teuchos::RCP<const NOX::Abstract::Vector> getNewtonPtr() const override { return dx_; }

 private:
  Teuchos::RCP<SolverFnBase<VectorClass> > fn_;
  Teuchos::RCP<NoxVector<VectorClass> > x_;
  Teuchos::RCP<NoxVector<VectorClass> > f_;
  Teuchos::RCP<NoxVector<VectorClass> > dx_;

  bool is_f_;
  bool is_jac_;
  bool is_dx_;
};

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
