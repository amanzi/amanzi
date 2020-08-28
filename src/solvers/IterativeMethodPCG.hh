/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Preconditioned conjugate gradient method for a linear solver.

/*!

.. _linear-solver-typed-pcg-spec:
.. admonition:: linear-solver-typed-pcg-spec

    * `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

    * `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

    * `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

    * `"convergence criterial`" ``[Array(string)]`` **"{relative rhs}"** A list of
      criteria, any of which can be applied.  Valid include:

      - `"relative rhs`" : measure error relative to the norm of the RHS vector
      - `"relative residual`" : measure error relative to the norm of the residual
      - `"absolute residual`" : measure error directly, norm of error
      - `"make one iteration`" : require at least one iteration to be performed before declaring success

*/

#ifndef AMANZI_PCG_OPERATOR_HH_
#define AMANZI_PCG_OPERATOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "InverseDefs.hh"
#include "InverseIterativeMethod.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix,
         class Preconditioner=Matrix,
         class Vector=typename Matrix::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
class IterativeMethodPCG :
      public InverseIterativeMethod<Matrix,Preconditioner,Vector,VectorSpace> {
 public:
  using InvIt = InverseIterativeMethod<Matrix,Preconditioner,Vector,VectorSpace>;

  IterativeMethodPCG() : InvIt() {}

  virtual int ApplyInverse(const Vector& v, Vector& hv) const override final {
    AMANZI_ASSERT(inited_ && h_.get());
    returned_code_ = PCG_(v, hv, tol_, max_itrs_, this->criteria_);
    if (returned_code_ <= 0) return 1;
    return 0;
  }

 protected:
  virtual std::string MethodName_() const override { return "PCG"; }

 private:
  int PCG_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

  IterativeMethodPCG(const IterativeMethodPCG& other) = delete;

 private:
  using InvIt::m_;
  using InvIt::h_;
  using InvIt::vo_;
  using InvIt::num_itrs_;
  using InvIt::max_itrs_;
  using InvIt::tol_;
  using InvIt::residual_;
  using InvIt::returned_code_;
  using InvIt::CheckConvergence_;
  using InvIt::inited_;
  using InvIt::rnorm0_;

};


/* ******************************************************************
 * PCG input/output data:
 *  f [input]         the right-hand side
 *  x [input/output]  initial guess / final solution
 *  tol [input]       convergence tolerance
 *  max_itrs [input]  maximum number of iterations
 *  criteria [input]  sum of termination critaria
 *
 *  Return value. If it is positive, it indicates the sucessful
 *  convergence criterion (criteria in a few exceptional cases) that
 *  was checked first. If it is negative, it indicates a failure, see
 *  LinearSolverDefs.hh for the error explanation.
 ***************************************************************** */
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
int IterativeMethodPCG<Matrix,Preconditioner,Vector,VectorSpace>::PCG_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Teuchos::OSTab tab = vo_->getOSTab();

  Vector r(f), p(f), v(f);  // construct empty vectors
  num_itrs_ = 0;

  double fnorm;
  f.Norm2(&fnorm);
  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    residual_ = 0.0;
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Converged, itr = " << num_itrs_ << " ||r|| = " << residual_ << std::endl;
    return criteria;  // Convergence for all criteria
  }

  m_->Apply(x, r);  // r = f - M * x
  r.Update(1.0, f, -1.0);
  r.Norm2(&rnorm0_);
  residual_ = rnorm0_;

  // Ignore all criteria if one iteration is enforced.
  if (! (criteria & LIN_SOLVER_MAKE_ONE_ITERATION)) {
    int ierr = CheckConvergence_(rnorm0_, fnorm);
    if (ierr != 0) return ierr;
  }

  // rare case that can happen in practise
  if (rnorm0_ == 0.0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Converged, itr = " << num_itrs_ << " ||r|| = " << residual_ << std::endl;
    return criteria;  // Convergence for all criteria
  }

  h_->ApplyInverse(r, p);  // gamma = (H r,r)
  double gamma0;
  p.Dot(r, &gamma0);
  if (gamma0 <= 0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Failed: non-SPD ApplyInverse: gamma0=" << gamma0 << std::endl;
    return LIN_SOLVER_NON_SPD_APPLY_INVERSE;
  }

  for (int i = 0; i < max_itrs; i++) {
    m_->Apply(p, v);
    double alpha;
    v.Dot(p, &alpha);

    if (alpha < 0.0) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        double pnorm;
        p.Norm2(&pnorm);
        *vo_->os() << "Failed: non-SPD Apply: alpha=" << alpha << " ||p||=" << pnorm << std::endl;
      }
      return LIN_SOLVER_NON_SPD_APPLY;
    }
    alpha = gamma0 / alpha;

    x.Update( alpha, p, 1.0);
    r.Update(-alpha, v, 1.0);

    h_->ApplyInverse(r, v);  // gamma1 = (H r, r)
    double gamma1;
    v.Dot(r, &gamma1);
    if (gamma1 < 0.0) {  // residual could be zero, so we use strict inequality
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Failed: non-SPD ApplyInverse: gamma1=" << gamma1 << std::endl;
      return LIN_SOLVER_NON_SPD_APPLY_INVERSE;
    }

    double rnorm;
    r.Norm2(&rnorm);
    residual_ = rnorm;

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << i << " ||r||=" << residual_ << std::endl;
    }
    num_itrs_ = i + 1;

    int ierr = CheckConvergence_(residual_, fnorm);
    if (ierr != 0) return ierr;

    double beta = gamma1 / gamma0;
    gamma0 = gamma1;

    p.Update(1.0, v, beta);
  }

  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "Failed (" << num_itrs_ << " itrs) ||r||=" << residual_ << " ||f||=" << fnorm << std::endl;
  return LIN_SOLVER_MAX_ITERATIONS;
};



}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
