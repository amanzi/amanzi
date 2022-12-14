/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Alicia Klinvex (amklinv@sandia.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! Trilinos/Belos implementations of iterative methods.
/*!
Includes GMRES

.. warning:: undocumented

*/

#pragma once
#include <cmath>

#include "BelosTypes.hpp"
#include "AmanziBelosOPWrapper.hh"
#include "AmanziBelosMVWrapper.hh"
#include "BelosPseudoBlockGmresSolMgr.hpp"

#include "Teuchos_RCP.hpp"
#include "errors.hh"
#include "VerboseObject.hh"

#include "InverseIterativeMethod.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Matrix,
          class Preconditioner = Matrix,
          class Vector = typename Matrix::Vector_t,
          class VectorSpace = typename Vector::Space_t>
class IterativeMethodBelos
  : public InverseIterativeMethod<Matrix, Preconditioner, Vector, VectorSpace> {
  using InvIt = InverseIterativeMethod<Matrix, Preconditioner, Vector, VectorSpace>;

 public:
  IterativeMethodBelos() : InvIt() {}

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final;
  virtual int ApplyInverse(const Vector& v, Vector& hv) const override final;

 protected:
  virtual std::string MethodName_() const override { return "Belos: GMRES"; }

 private:
  using InvIt::m_;
  using InvIt::h_;
  using InvIt::vo_;
  using InvIt::num_itrs_;
  using InvIt::max_itrs_;
  using InvIt::tol_;
  using InvIt::criteria_;
  using InvIt::residual_;
  using InvIt::returned_code_;
  using InvIt::krylov_dim_;
  using InvIt::inited_;

  Teuchos::RCP<Teuchos::ParameterList> belos_list_;
};


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template <class Matrix, class Preconditioner, class Vector, class VectorSpace>
void
IterativeMethodBelos<Matrix, Preconditioner, Vector, VectorSpace>::set_inverse_parameters(
  Teuchos::ParameterList& plist)
{
  InvIt::set_inverse_parameters(plist);

  belos_list_ = Teuchos::rcp(new Teuchos::ParameterList());
  belos_list_->set("Num Blocks", krylov_dim_);
  belos_list_->set("Maximum Iterations", max_itrs_);
  belos_list_->set("Maximum Restarts", 2 * krylov_dim_ * max_itrs_);
  belos_list_->set("Convergence Tolerance", tol_);
  belos_list_->set("Output Frequency", 1);

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    belos_list_->set("Verbosity",
                     Belos::Errors + Belos::Warnings + Belos::IterationDetails +
                       Belos::StatusTestDetails + Belos::Debug + Belos::TimingDetails);
  } else if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    belos_list_->set("Verbosity",
                     Belos::Errors + Belos::Warnings + Belos::IterationDetails +
                       Belos::StatusTestDetails);
  } else if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    belos_list_->set("Verbosity", Belos::Errors + Belos::Warnings);
  } else if (vo_->os_OK(Teuchos::VERB_LOW)) {
    belos_list_->set("Verbosity", Belos::Errors + Belos::Warnings);
  }

  if (criteria_ & LIN_SOLVER_RELATIVE_RHS) {
    belos_list_->set("Implicit Residual Scaling", "Norm of RHS");
    belos_list_->set("Explicit Residual Scaling", "Norm of RHS");
  } else if (criteria_ & LIN_SOLVER_RELATIVE_RESIDUAL) {
    belos_list_->set("Implicit Residual Scaling", "Norm of Initial Residual");
    belos_list_->set("Explicit Residual Scaling", "Norm of Initial Residual");
  } else {
    belos_list_->set("Implicit Residual Scaling", "None");
    belos_list_->set("Explicit Residual Scaling", "None");
  }
}

template <class Matrix, class Preconditioner, class Vector, class VectorSpace>
int
IterativeMethodBelos<Matrix, Preconditioner, Vector, VectorSpace>::ApplyInverse(const Vector& v,
                                                                                Vector& hv) const
{
  typedef Belos::LinearProblem<double, Belos::MultiVec<double>, Belos::Operator<double>>
    LinearProblem;
  typedef Belos::PseudoBlockGmresSolMgr<double, Belos::MultiVec<double>, Belos::Operator<double>>
    GmresSolver;

  AMANZI_ASSERT(belos_list_.get()); // set_inverse_parameters() called
  AMANZI_ASSERT(m_.get());          // set_matrices() called

  // NOTE: this is not efficiently implemented and should get fixed.  The
  // problem can get created and reused, but not in the current pattern. It may
  // be as simple as keeping the LinearProblem from call to call.

  auto mat = Teuchos::rcp(new AmanziBelosOp<Matrix, Vector>(this->m_));
  auto prec = Teuchos::rcp(new AmanziBelosOp<Preconditioner, Vector>(this->h_, true));
  auto lhs = Teuchos::rcp(new CompositeMultiVector<Vector>(Teuchos::rcp(&hv, false)));
  auto rhs =
    Teuchos::rcp(new CompositeMultiVector<Vector>(Teuchos::rcp(const_cast<Vector*>(&v), false)));
  auto problem = Teuchos::rcp(new LinearProblem(mat, lhs, rhs));

  problem->setLeftPrec(prec);
  problem->setProblem();

  GmresSolver solver(problem, belos_list_);
  Belos::ReturnType success = solver.solve();
  num_itrs_ = solver.getNumIters();
  residual_ = solver.achievedTol();

  if (success == Belos::Converged) {
    returned_code_ = LIN_SOLVER_BELOS_SAYS_SUCCESS;
    return 0;
  } else {
    returned_code_ = LIN_SOLVER_BELOS_SAYS_FAIL;
    return 1;
  }
}

} // namespace AmanziSolvers
} // namespace Amanzi
