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

Parameters:

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

#include "VerboseObject.hh"

#include "LinearOperatorDefs.hh"
#include "LinearOperator.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorPCG : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorPCG(const Teuchos::RCP<const Matrix>& m,
                    const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h),
      tol_(1e-6),
      overflow_tol_(3.0e+50),  // mass of the Universe (J.Hopkins)
      max_itrs_(100),
      criteria_(LIN_SOLVER_RELATIVE_RHS),
      initialized_(false) {};

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix,Vector,VectorSpace>::Init(); }

  int ApplyInverse(const Vector& v, Vector& hv) const {
    if (!initialized_) {
      Errors::Message msg("LinearOperatorPCG: has not been initialized.");
      Exceptions::amanzi_throw(msg);
    }
    int ierr = PCG_(v, hv, tol_, max_itrs_, criteria_);
    returned_code_ = ierr;
    //    return (ierr > 0) ? 0 : 1;  // Positive ierr code means success.
    return returned_code_;
  }

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void add_criteria(int criteria) { criteria_ |= criteria; }
  void set_overflow(double tol) { overflow_tol_ = tol; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }
  int returned_code() { return returned_code_; }

 public:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  int PCG_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

  LinearOperatorPCG(const LinearOperatorPCG& other); // not implemented

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::h_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  int max_itrs_, criteria_;
  double tol_, overflow_tol_;
  mutable int num_itrs_, returned_code_;
  mutable double residual_;
  mutable bool initialized_;
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
template<class Matrix, class Vector, class VectorSpace>
int LinearOperatorPCG<Matrix, Vector, VectorSpace>::PCG_(
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
  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;

  // Ignore all criteria if one iteration is enforced.
  if (rnorm0 > overflow_tol_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Diverged, ||r||=" << rnorm0 << std::endl;
    return LIN_SOLVER_RESIDUAL_OVERFLOW;
  }

  if (! (criteria & LIN_SOLVER_MAKE_ONE_ITERATION)) {
    if (criteria & LIN_SOLVER_RELATIVE_RHS) {
      if (rnorm0 < tol * fnorm) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative RHS), itr=" << num_itrs_ 
                     << " ||r||=" << rnorm0 << " ||f||=" << fnorm << std::endl;
        return LIN_SOLVER_RELATIVE_RHS;
      }
    } else if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm0 < tol) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (absolute res), itr=" << num_itrs_ << "||r||=" << rnorm0 << std::endl;
        return LIN_SOLVER_ABSOLUTE_RESIDUAL;
      }
    }
  }

  // rare case that can happen in practise
  if (rnorm0 == 0.0) {
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

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << i << " ||r||=" << residual_ << std::endl;
    }

    if (rnorm > overflow_tol_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Diverged, ||r||=" << rnorm << std::endl;
      return LIN_SOLVER_RESIDUAL_OVERFLOW;
    }

    // Return the first criterion which is fulfilled.
    num_itrs_ = i + 1;
    if (criteria & LIN_SOLVER_RELATIVE_RHS) {
      if (rnorm < tol * fnorm) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative RHS), itr=" << num_itrs_ 
                     << " ||r||=" << rnorm << " ||f||=" << fnorm << std::endl;
        return LIN_SOLVER_RELATIVE_RHS;
      }
    } else if (criteria & LIN_SOLVER_RELATIVE_RESIDUAL) {
      if (rnorm < tol * rnorm0) { 
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative res), itr=" << num_itrs_ 
                     << " ||r||=" << rnorm << " ||r0||=" << rnorm0 << std::endl;
        return LIN_SOLVER_RELATIVE_RESIDUAL;
      }
    } else if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm < tol) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (absolute res), itr=" << num_itrs_ << " ||r||=" << rnorm << std::endl;
        return LIN_SOLVER_ABSOLUTE_RESIDUAL;
      }
    }

    double beta = gamma1 / gamma0;
    gamma0 = gamma1;

    p.Update(1.0, v, beta);
  }

  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "Failed (" << num_itrs_ << " itrs) ||r||=" << residual_ << " ||f||=" << fnorm << std::endl;

  return LIN_SOLVER_MAX_ITERATIONS;
};


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorPCG<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::PCG", plist));

  tol_ = plist.get<double>("error tolerance", 1e-6);
  max_itrs_ = plist.get<int>("maximum number of iterations", 100);
  overflow_tol_ = plist.get<double>("overflow tolerance", 3.0e+50);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names = plist.get<Teuchos::Array<std::string> > ("convergence criteria").toVector();

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
	msg << "LinearOperatorGMRES: \"convergence criteria\" type \"" << names[i] << "\" is not recognized.";
	Exceptions::amanzi_throw(msg);
      }
    }
  } else {
    criteria = LIN_SOLVER_RELATIVE_RHS;
  }

  set_criteria(criteria);
  initialized_ = true;
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
