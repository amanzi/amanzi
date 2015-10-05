/*
  num_itrs_accumulated_ = 0;
  This is the Linear Solver component of the Amanzi code.
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
  Konstantin Lipnikov (lipnikov@lanl.gov)

  Generalized minimum residual method (Yu.Kuznetsov, 1968; Y.Saad, 1986)
  Usage:
*/

#ifndef  AMANZI_GMRES_OPERATOR_HH_
#define  AMANZI_GMRES_OPERATOR_HH_

#include <cmath>

#include "Teuchos_RCP.hpp"
#include "errors.hh"
#include "VerboseObject.hh"

#include "DenseMatrix.hh"
#include "LinearOperator.hh"
#include "LinearOperatorDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorGMRES : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorGMRES(const Teuchos::RCP<const Matrix>& m,
                      const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h),
      tol_(1e-6),
      overflow_tol_(3.0e+50),  // mass of the Universe (J.Hopkins)
      max_itrs_(100),
      krylov_dim_(10),
      criteria_(LIN_SOLVER_RELATIVE_RHS),
      initialized_(false) {}

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix, Vector, VectorSpace>::Init(); }

  int ApplyInverse(const Vector& v, Vector& hv) const {
    if (!initialized_) {
      Errors::Message msg("LinearOperatorGMRES: has not been initialized.");
      Exceptions::amanzi_throw(msg);
    }

    int ierr = GMRESRestart_(v, hv, tol_, max_itrs_, criteria_);
    returned_code_ = ierr;
    //return (ierr > 0) ? 0 : 1;  // Positive ierr code means success.
    return returned_code_;
  }

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void add_criteria(int criteria) { criteria_ |= criteria; }
  void set_krylov_dim(int n) { krylov_dim_ = n; }
  void set_overflow(double tol) { overflow_tol_ = tol; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }
  int returned_code() { return returned_code_; }

 public:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  int GMRESRestart_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  int GMRES_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  void ComputeSolution_(Vector& x, int k, WhetStone::DenseMatrix& T, double* s,
                        std::vector<Teuchos::RCP<Vector> >& v) const;
  void InitGivensRotation_( double& dx, double& dy, double& cs, double& sn) const;
  void ApplyGivensRotation_(double& dx, double& dy, double& cs, double& sn) const;

  LinearOperatorGMRES(const LinearOperatorGMRES& other); // not implemented

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::h_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  int max_itrs_, criteria_, krylov_dim_;
  double tol_, overflow_tol_;
  mutable int num_itrs_, num_itrs_total_, returned_code_;
  mutable double residual_, fnorm_, rnorm0_;
  mutable bool initialized_;

  int controller_start_, controller_end_;
  mutable double controller_[2];
};


/* ******************************************************************
 * GMRES with restart input/output data:
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
int LinearOperatorGMRES<Matrix, Vector, VectorSpace>::GMRESRestart_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  num_itrs_total_ = 0;
  rnorm0_ = -1.0;

  int ierr(LIN_SOLVER_MAX_ITERATIONS);
  while (ierr == LIN_SOLVER_MAX_ITERATIONS && num_itrs_total_ < max_itrs) {
    int max_itrs_left = max_itrs - num_itrs_total_;

    ierr = GMRES_(f, x, tol, max_itrs_left, criteria);
    if (ierr == LIN_SOLVER_RESIDUAL_OVERFLOW) return ierr;
  }

  if (ierr == LIN_SOLVER_MAX_ITERATIONS) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Not converged (max iterations), ||r||=" << residual_ 
                 << " ||f||=" << fnorm_ << std::endl;
  }

  num_itrs_ = num_itrs_total_;
  return ierr;
}


/* ******************************************************************
 * GMRES input/output data:
 *  f [input]         the right-hand side
 *  x [input/output]  initial guess / final solution
 *  tol [input]       convergence tolerance
 *  max_itrs [input]  maximum number of iterations
 *  criteria [input]  sum of termination critaria
 *
 *  Return value. See above.
 ***************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
int LinearOperatorGMRES<Matrix, Vector, VectorSpace>::GMRES_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Teuchos::OSTab tab = vo_->getOSTab();

  Vector w(f), r(f), p(f);  // construct empty vectors
  std::vector<Teuchos::RCP<Vector> > v(krylov_dim_ + 1);

  double s[krylov_dim_ + 1], cs[krylov_dim_ + 1], sn[krylov_dim_ + 1];
  WhetStone::DenseMatrix T(krylov_dim_ + 1, krylov_dim_);
  num_itrs_ = 0;

  double fnorm;
  // h_->ApplyInverse(f, r);
  // r.Dot(fnorm, f);  This is the preconditioned norm of the residual.
  f.Norm2(&fnorm);
  fnorm_ = fnorm;

  m_->Apply(x, p);  // p = f - M * x
  p.Update(1.0, f, -1.0);

  h_->ApplyInverse(p, r);
  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;
  if (rnorm0_ < 0.0) rnorm0_ = rnorm0;

  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Converged, itr=" << num_itrs_total_ << " ||r||=" << rnorm0 << std::endl;
    return criteria;  // Zero solution satifies all criteria.
  }

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
	  *vo_->os() << "Converged (relative RHS), itr=" << num_itrs_total_ 
                     << " ||r||=" << rnorm0 << std::endl;
	return LIN_SOLVER_RELATIVE_RHS;
      }
    } else if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm0 < tol) {
	if (vo_->os_OK(Teuchos::VERB_MEDIUM))
	  *vo_->os() << "Converged (absolute res), itr=" << num_itrs_total_ 
                     << " ||r||=" << rnorm0 << std::endl;
	return LIN_SOLVER_ABSOLUTE_RESIDUAL;
      }
    }
  }

  v[0] = Teuchos::rcp(new Vector(r));
  v[0]->Update(0.0, r, 1.0 / rnorm0);

  s[0] = rnorm0;

  for (int i = 0; i < krylov_dim_; i++) {
    m_->Apply(*(v[i]), p);
    h_->ApplyInverse(p, w);

    double tmp(0.0);
    for (int k = 0; k <= i; k++) {  // Arnoldi algorithm
      w.Dot(*(v[k]), &tmp);
      w.Update(-tmp, *(v[k]), 1.0);
      T(k, i) = tmp;
    }
    w.Norm2(&tmp);
    T(i + 1, i) = tmp;
    s[i + 1] = 0.0;

    for (int k = 0; k < i; k++) {
      ApplyGivensRotation_(T(k, i), T(k + 1, i), cs[k], sn[k]);
    }

    InitGivensRotation_( T(i, i), T(i + 1, i), cs[i], sn[i]);
    ApplyGivensRotation_(T(i, i), T(i + 1, i), cs[i], sn[i]);
    ApplyGivensRotation_(s[i],    s[i + 1],    cs[i], sn[i]);
    residual_ = fabs(s[i + 1]);

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << num_itrs_total_ << " ||r||=" << residual_ << std::endl;
    }

    // Check all criteria one-by-one.
    num_itrs_ = i + 1;
    num_itrs_total_++;
    if (criteria & LIN_SOLVER_RELATIVE_RHS) {
      if (residual_ < tol * fnorm) {
        ComputeSolution_(x, i, T, s, v);  // vector s is overwritten
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative RHS), itr=" << num_itrs_total_ 
                     << " ||r||=" << residual_ << " ||f||=" << fnorm << std::endl;
        return LIN_SOLVER_RELATIVE_RHS;
      }
    }

    if (criteria & LIN_SOLVER_RELATIVE_RESIDUAL) {
      if (residual_ < tol * rnorm0_) {
        ComputeSolution_(x, i, T, s, v);  // vector s is overwritten
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (relative res), itr=" << num_itrs_total_ 
                     << " ||r||=" << residual_ << " ||r0||=" << rnorm0_ << std::endl;
        return LIN_SOLVER_RELATIVE_RESIDUAL;
      }
    }

    if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (residual_ < tol) {
        ComputeSolution_(x, i, T, s, v);  // vector s is overwritten
        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
          *vo_->os() << "Converged (absolute res), itr=" << num_itrs_total_ 
                     << " ||r||=" << residual_ << std::endl;
        return LIN_SOLVER_ABSOLUTE_RESIDUAL;
      }
    }

    // optional controller of convergence
    if (i == controller_start_) { 
      controller_[0] = residual_;
    } else if (i == controller_end_) {
      double len = 0.5 / (controller_end_ - controller_start_);
      controller_[0] = std::pow(controller_[0] / residual_, len);
      controller_[0] = std::min(controller_[0], 2.0);
      controller_[1] = residual_;
    } else if (i > controller_end_) {
      double reduction = controller_[1] / residual_;
      if (reduction < controller_[0]) {
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "controller indicates convergence stagnation\n"; 

        ComputeSolution_(x, i, T, s, v); 
        return LIN_SOLVER_MAX_ITERATIONS;
      }
      controller_[1] = residual_;
    }

    if (i < krylov_dim_ - 1) {
      v[i + 1] = Teuchos::rcp(new Vector(w));
      if (tmp != 0.0) {  // zero occurs in exact arithmetic
        v[i + 1]->Update(0.0, r, 1.0 / tmp);
      }
    }
  }

  ComputeSolution_(x, krylov_dim_ - 1, T, s, v);  // vector s is overwritten
  return LIN_SOLVER_MAX_ITERATIONS;
}


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorGMRES<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::GMRES", plist));

  tol_ = plist.get<double>("error tolerance", 1e-16);
  max_itrs_ = plist.get<int>("maximum number of iterations", 100);
  krylov_dim_ = plist.get<int>("size of Krylov space", 10);
  overflow_tol_ = plist.get<double>("overflow tolerance", 3.0e+50);

  controller_start_ = plist.get<int>("controller training start", 0);
  controller_end_ = plist.get<int>("controller training end", 3);
  controller_end_ = std::max(controller_end_, controller_start_ + 1);

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
      }
    }
  } else {
    criteria = LIN_SOLVER_RELATIVE_RHS;
  }

  set_criteria(criteria);
  initialized_ = true;
}


/* ******************************************************************
 * Givens rotations: initialization
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorGMRES<Matrix, Vector, VectorSpace>::InitGivensRotation_(
    double& dx, double& dy, double& cs, double& sn) const
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    double tmp = dx / dy;
    sn = 1.0 / sqrt(1.0 + tmp * tmp);
    cs = tmp * sn;
  } else {
    double tmp = dy / dx;
    cs = 1.0 / sqrt(1.0 + tmp * tmp);
    sn = tmp * cs;
  }
}


/* ******************************************************************
 * Givens rotations: applications
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorGMRES<Matrix, Vector, VectorSpace>::ApplyGivensRotation_(
    double& dx, double& dy, double& cs, double& sn) const
{
  double tmp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = tmp;
}


/* ******************************************************************
 * Computation of the solution destroys vector s.
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorGMRES<Matrix, Vector, VectorSpace>::ComputeSolution_(
    Vector& x, int k, WhetStone::DenseMatrix& T, double* s,
    std::vector<Teuchos::RCP<Vector> >& v) const
{
  for (int i = k; i >= 0; i--) {
    s[i] /= T(i, i);

    for (int j = i - 1; j >= 0; j--) {
      s[j] -= T(j, i) * s[i];
    }
  }

  for (int j = 0; j <= k; j++) {
    x.Update(s[j], *(v[j]), 1.0);
  }
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
