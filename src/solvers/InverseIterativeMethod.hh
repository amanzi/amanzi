/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Base class for iterative methods, including Krylov methods.
/*!

Document me...

*/
#pragma once

#include "Key.hh"
#include "VerboseObject.hh"
#include "Inverse.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Matrix,
          class Preconditioner = Matrix,
          class Vector = typename Matrix::Vector_t,
          class VectorSpace = typename Vector::Space_t>
class InverseIterativeMethod : public Inverse<Matrix, Preconditioner, Vector, VectorSpace> {
 public:
  InverseIterativeMethod()
    : Inverse<Matrix, Preconditioner, Vector, VectorSpace>(), inited_(false){};

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override;
  virtual void InitializeInverse() override { h_->InitializeInverse(); }
  virtual void ComputeInverse() override { h_->ComputeInverse(); }

  // control and statistics -- must be valid for both iterative and
  // non-iterative methods, approximate and exact methods.
  virtual double residual() const override { return residual_; }
  virtual int num_itrs() const override { return num_itrs_; }

  virtual void add_criteria(int criteria) override { criteria_ |= criteria; }
  void set_criteria(int criteria)
  {
    criteria_ = criteria;
    if (!(criteria & (LIN_SOLVER_RELATIVE_RHS | LIN_SOLVER_RELATIVE_RESIDUAL |
                      LIN_SOLVER_ABSOLUTE_RESIDUAL))) {
      // need at least one of these three
      Errors::Message msg("InverseIterativeMethod: criteria must include one of RELATIVE_RHS, "
                          "RELATIVE_RESIDUAL, or ABSOLUTE_RESIDUAL");
      Exceptions::amanzi_throw(msg);
    }
  }

  virtual int returned_code() const override { return returned_code_; }

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_krylov_dim(int n) { krylov_dim_ = n; }
  void set_overflow(double tol) { overflow_tol_ = tol; }

  virtual std::string returned_code_string() const override
  {
    switch (this->returned_code()) {
    case AmanziSolvers::LIN_SOLVER_NON_SPD_APPLY:
      return "Linear system is not SPD.";
    case AmanziSolvers::LIN_SOLVER_NON_SPD_APPLY_INVERSE:
      return "Linear system is not SPD.";
    case AmanziSolvers::LIN_SOLVER_MAX_ITERATIONS:
      return "Maximum iterations are reached in solution of linear system.";
    case AmanziSolvers::LIN_SOLVER_RESIDUAL_OVERFLOW:
      return "Residual overflow in solution of linear system.";
    }

    if (this->returned_code() < 0) {
      return "Iterative method returned an unrecoverable error code: " +
             std::to_string(this->returned_code()) + ".";
    } else if (this->returned_code() == 0) {
      return "Iterative method is still iterating.";
    } else {
      return "Iterative method converged in " + std::to_string(this->returned_code()) +
             " iterations.";
    }
  }

 protected:
  int CheckConvergence_(double rnorm, double fnorm) const;
  virtual std::string MethodName_() const = 0;

 protected:
  int max_itrs_, criteria_, krylov_dim_;
  double tol_, overflow_tol_;
  mutable int num_itrs_, returned_code_;
  mutable double residual_, rnorm0_;
  bool inited_;

  using Inverse<Matrix, Preconditioner, Vector, VectorSpace>::m_;
  using Inverse<Matrix, Preconditioner, Vector, VectorSpace>::h_;
  Teuchos::RCP<VerboseObject> vo_;
};


//
// Parse the input spec.
//
template <class Matrix, class Preconditioner, class Vector, class VectorSpace>
void inline InverseIterativeMethod<Matrix, Preconditioner, Vector, VectorSpace>::
  set_inverse_parameters(Teuchos::ParameterList& plist)
{
  std::string vo_name = this->name() + "::" + this->MethodName_();
  vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist));

  tol_ = plist.template get<double>("error tolerance", 1e-16);
  max_itrs_ = plist.template get<int>("maximum number of iterations", 100);
  krylov_dim_ = plist.template get<int>("size of Krylov space", 10);
  overflow_tol_ = plist.template get<double>("overflow tolerance", 3.0e+50);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names = plist.template get<Teuchos::Array<std::string>>("convergence criteria").toVector();

    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "relative rhs") {
        criteria += LIN_SOLVER_RELATIVE_RHS;
      } else if (names[i] == "relative residual") {
        criteria += LIN_SOLVER_RELATIVE_RESIDUAL;
      } else if (names[i] == "absolute residual") {
        criteria += LIN_SOLVER_ABSOLUTE_RESIDUAL;
      } else if (names[i] == "make one iteration") {
        criteria += LIN_SOLVER_MAKE_ONE_ITERATION;
      } else {
        Errors::Message msg;
        msg << "InverseIterativeMethod: \"convergence criteria\" type \"" << names[i]
            << "\" is not recognized.";
        Exceptions::amanzi_throw(msg);
      }
    }
  } else {
    criteria = LIN_SOLVER_RELATIVE_RHS;
  }

  set_criteria(criteria);
  inited_ = true;
}


//
// Check whether convergence criteria are satisfied.
//
template <class Matrix, class Preconditioner, class Vector, class VectorSpace>
int inline InverseIterativeMethod<Matrix, Preconditioner, Vector, VectorSpace>::CheckConvergence_(
  double rnorm,
  double fnorm) const
{
  if (rnorm > overflow_tol_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) *vo_->os() << "Diverged, ||r||=" << rnorm << std::endl;
    return LIN_SOLVER_RESIDUAL_OVERFLOW;
  }

  if (criteria_ & LIN_SOLVER_RELATIVE_RHS) {
    if (rnorm < tol_ * fnorm) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Converged (relative RHS), itr=" << num_itrs_ << " ||r||=" << rnorm
                   << " ||f||=" << fnorm << std::endl;
      return LIN_SOLVER_RELATIVE_RHS;
    }
  }

  if (criteria_ & LIN_SOLVER_RELATIVE_RESIDUAL) {
    if (rnorm < tol_ * rnorm0_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Converged (relative res), itr=" << num_itrs_ << " ||r||=" << residual_
                   << " ||r0||=" << rnorm0_ << std::endl;
      return LIN_SOLVER_RELATIVE_RESIDUAL;
    }
  }

  if (criteria_ & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
    if (rnorm < tol_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Converged (absolute res), itr=" << num_itrs_ << " ||r||=" << rnorm
                   << std::endl;
      return LIN_SOLVER_ABSOLUTE_RESIDUAL;
    }
  }

  return 0;
};


} // namespace AmanziSolvers
} // namespace Amanzi
