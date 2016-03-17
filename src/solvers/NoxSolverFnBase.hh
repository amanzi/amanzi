#ifndef AMANZI_NOX_SOLVER_FN_BASE_
#define AMANZI_NOX_SOLVER_FN_BASE_

#include "NOX_Abstract_Group.H"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector>
class NoxSolverFnBase : public NOX::Abstract::Group {

 public:
  NoxSolverFnBase(const Teuchos::RCP<SolverFnBase<Vector> >& fn) :
      fn_(fn) {}

  NoxSolverFnBase<Vector>&
  operator=(const NoxSolverFnBase<Vector>& source) {
    ASSERT(0);
  }

  NOX::Abstract::Group::ReturnType
  computeF() {
    fn_->Residual(x_->get_vector(), f_->get_vector());
    return NOX::Abstract::Group::Ok;
  }

  NOX::Abstract::Group::ReturnType
  computeJacobian() {
    fn_->UpdatePreconditioner(x_->get_vector());
    return NOX::Abstract::Group::Ok;
  }

  // NOTE: Room for improvement for things like Eisenstat-Walker for tuning
  // tols by passing these params into a PK's linear operator
  NOX::Abstract::Group::ReturnType
  computeNewton(Teuchos::ParameterList& params) {
    ASSERT(isJacobian());
    
    int ierr = fn_->ApplyPreconditioner(f_->get_vector(), dx_->get_vector());
    dx_->Scale(-1.);
    if (ierr > 0) {
      return NOX::Abstract::Group::Ok;
    } else {
      return NOX::Abstract::Group::NotConverged;
    }
  }

  
  NOX::Abstract::Group::ReturnType
  applyJacobianInverse(Teuchos::ParameterList& params,
                       const NOX::Abstract::Vector& input,
                       NOX::Abstract::Vector& result) {
    const NoxVector<Vector>* input_a =
        std::dynamic_cast<NoxVector<Vector>*>(&input);
    const NoxVector<Vector>* result_a =
        std::dynamic_cast<NoxVector<Vector>*>(&result);

    int ierr = fn_->ApplyPreconditioner(input_a.get_vector(), result_a.get_vector());
    dx_->Scale(-1.);
    if (ierr > 0) {
      return NOX::Abstract::Group::Ok;
    } else {
      return NOX::Abstract::Group::NotConverged;
    }
  }

  void setX(x) {
    fn_->ChangedSolution();
    x_ = ..
  }
  
 private:
  Teuchos::RCP<SolverFnBase<Vector> > fn_;
  Teuchos::RCP<NoxVector<Vector> > x_;
  Teuchos::RCP<NoxVector<Vector> > f_;
  Teuchos::RCP<NoxVector<Vector> > dx_;

  bool is_f_;
  bool is_jac_;
  bool is_dx_;
  
};



} // namespace
} // namespace


#endif
