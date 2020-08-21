/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Generalized minimum residual method for a linear solver.

/*!

Based on the methods of Yu. Kuznetsov, 1968; Y.Saad, 1986.  Deflated version of
GMRES is due to R.Morgan, GMRES with deflated restarting, 2002 SISC; S.Rollin,
W.Fichtner, Improving accuracy of GMRES with deflated restarting, 2007 SISC.

.. _linear-solver-typed-gmres-spec:
.. admonition:: linear-solver-typed-gmres-spec

    * `"error tolerance`" ``[double]`` **1.e-6** Tolerance on which to declare success.

    * `"maximum number of iterations`" ``[int]`` **100** Maximum iterations before declaring failure.

    * `"overflow tolerance`" ``[double]`` **3.e50** Error above this value results in failure.

    * `"convergence criterial`" ``[Array(string)]`` **"{relative rhs}"** A list of
      criteria, any of which can be applied.  Valid include:

      - `"relative rhs`" : measure error relative to the norm of the RHS vector
      - `"relative residual`" : measure error relative to the norm of the residual
      - `"absolute residual`" : measure error directly, norm of error
      - `"make one iteration`" : require at least one iteration to be performed before declaring success

    * `"size of Krylov space`" ``[int]`` **10** Size of the Krylov space used to span the residual.

    * `"controller training start`" ``[int]`` **0** Start iteration for determining
      convergence rates. (Add more please!)

    * `"controller training end`" ``[int]`` **3** Start iteration for determining
      convergence rates. (Add more please!)

    * `"preconditioning strategy`" ``[string]`` **left** Valid are "left" and
      "right"-type preconditioning (see Saad 1986)

    * `"maximum size of deflation space`" ``[int]`` **0** Size of the deflation space, see Rollin et al.

*/

#ifndef AMANZI_GMRES_OPERATOR_HH_
#define AMANZI_GMRES_OPERATOR_HH_

#include <cmath>

#include "Teuchos_RCP.hpp"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "InverseIterativeMethod.hh"
#include "InverseDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix,
         class Preconditioner=Matrix,
         class Vector=typename Matrix::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
class IterativeMethodGMRES :
      public InverseIterativeMethod<Matrix,Preconditioner,Vector,VectorSpace> {
 private:
  using InvIt = InverseIterativeMethod<Matrix,Preconditioner,Vector,VectorSpace>;

 public:
  IterativeMethodGMRES() :
      InvIt() {}

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final;

  virtual int ApplyInverse(const Vector& v, Vector& hv) const override final {
    AMANZI_ASSERT(inited_ && h_.get());
    returned_code_ = GMRESRestart_(v, hv, this->tol_, this->max_itrs_, this->criteria_);
    return returned_code_;
  }
  
 protected:
  virtual std::string MethodName_() const override { return "GMRES"; }

 private:
  int GMRESRestart_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  int GMRES_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  int GMRES_Deflated_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

  void ComputeSolution_(Vector& x, int k, WhetStone::DenseMatrix& T, double* s,
                        Vector& p, Vector& r) const;
  void ComputeSolution_(Vector& x, double* d, Vector& p, Vector& r) const;

  void InitGivensRotation_( double& dx, double& dy, double& cs, double& sn) const;
  void ApplyGivensRotation_(double& dx, double& dy, double& cs, double& sn) const;
 
  IterativeMethodGMRES(const IterativeMethodGMRES& other) = delete;

 private:
  using InvIt::m_;
  using InvIt::h_;
  using InvIt::vo_;
  using InvIt::num_itrs_;
  using InvIt::residual_;
  using InvIt::returned_code_;
  using InvIt::CheckConvergence_;
  using InvIt::krylov_dim_;
  using InvIt::inited_;

  mutable std::vector<Teuchos::RCP<Vector> > v_;
  mutable WhetStone::DenseMatrix Hu_;  // upper Hessenberg matrix

  mutable double rnorm0_, fnorm_;
  
  mutable int num_itrs_inner_;
  int controller_start_, controller_end_;
  mutable double controller_[2];

  bool left_pc_;
  int deflation_;
  mutable int num_ritz_;
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
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
int IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::GMRESRestart_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  // initialize verbose object
  Teuchos::OSTab tab = vo_->getOSTab();

  // allocate memory for Krylov space
  v_.resize(krylov_dim_ + 1, Teuchos::null);

  num_itrs_ = 0;
  rnorm0_ = -1.0;

  int ierr(LIN_SOLVER_MAX_ITERATIONS);
  while (ierr == LIN_SOLVER_MAX_ITERATIONS && num_itrs_ < max_itrs) {
    int max_itrs_left = max_itrs - num_itrs_;

    if (deflation_ == 0) {
      ierr = GMRES_(f, x, tol, max_itrs_left, criteria);
    } else {
      ierr = GMRES_Deflated_(f, x, tol, max_itrs_left, criteria);
    }
    if (ierr == LIN_SOLVER_RESIDUAL_OVERFLOW) return ierr;
  }

  if (ierr == LIN_SOLVER_MAX_ITERATIONS) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Not converged (max iterations), ||r||=" << residual_ 
                 << " ||f||=" << fnorm_ << std::endl;
  }

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
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
int IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::GMRES_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Vector w(f), r(f), p(f);  // construct empty vectors

  double s[krylov_dim_ + 1], cs[krylov_dim_ + 1], sn[krylov_dim_ + 1];
  WhetStone::DenseMatrix T(krylov_dim_ + 1, krylov_dim_);
  num_itrs_inner_ = 0;

  double fnorm;
  // h_->ApplyInverse(f, r);
  // r.Dot(fnorm, f);  This is the preconditioned norm of the residual.
  f.Norm2(&fnorm);
  fnorm_ = fnorm;

  // initial residual is r = f - M x for the right preconditioner
  // and r = H (f - M x)  for the left preconditioner
  if (left_pc_) {
    m_->Apply(x, p);
    p.Update(1.0, f, -1.0);
    h_->ApplyInverse(p, r);
  } else {
    m_->Apply(x, r);
    r.Update(1.0, f, -1.0);
  }

  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;
  if (rnorm0_ < 0.0) rnorm0_ = rnorm0;

  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Converged, itr=" << num_itrs_ << " ||r||=" << rnorm0 << std::endl;
    return criteria;  // Zero solution satifies all criteria.
  }

  // Ignore all criteria if one iteration is enforced.
  if (! (criteria & LIN_SOLVER_MAKE_ONE_ITERATION)) {
    int ierr = CheckConvergence_(rnorm0, fnorm);
    if (ierr != 0) return ierr;
  }

  v_[0] = Teuchos::rcp(new Vector(r));
  v_[0]->Update(0.0, r, 1.0 / rnorm0);

  s[0] = rnorm0;

  for (int i = 0; i < krylov_dim_; i++) {
    // calculate H M v_i for the left preconditioner
    //       and M H v_i for the right preconditioner
    if (left_pc_) {
      m_->Apply(*(v_[i]), p);
      h_->ApplyInverse(p, w);
    } else {
      h_->ApplyInverse(*(v_[i]), p);
      m_->Apply(p, w);
    }

    double tmp(0.0);
    for (int k = 0; k <= i; k++) {  // Arnoldi algorithm
      w.Dot(*(v_[k]), &tmp);
      w.Update(-tmp, *(v_[k]), 1.0);
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

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << num_itrs_ << " ||r||=" << residual_ << std::endl;
    }

    // Check all criteria one-by-one.
    num_itrs_inner_ = i + 1;
    num_itrs_++;

    int ierr = CheckConvergence_(residual_, fnorm);
    if (ierr != 0) {
      ComputeSolution_(x, i, T, s, p, r);  // vector s is overwritten
      return ierr;
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

        ComputeSolution_(x, i, T, s, p, r); 
        return LIN_SOLVER_MAX_ITERATIONS;
      }
      controller_[1] = residual_;
    }

    if (i < krylov_dim_ - 1) {
      v_[i + 1] = Teuchos::rcp(new Vector(w));
      if (tmp != 0.0) {  // zero occurs in exact arithmetic
        v_[i + 1]->Update(0.0, r, 1.0 / tmp);
      }
    }
  }

  ComputeSolution_(x, krylov_dim_ - 1, T, s, p, r);  // vector s is overwritten
  return LIN_SOLVER_MAX_ITERATIONS;
}


/* ******************************************************************
* GMRES with deflated start, input/output data: see GMRES_(...)
***************************************************************** */
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
int IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::GMRES_Deflated_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Vector p(f), r(f), w(f);
  WhetStone::DenseVector d(krylov_dim_ + 1), g(krylov_dim_);
  WhetStone::DenseMatrix T(krylov_dim_ + 1, krylov_dim_);

  double fnorm;
  f.Norm2(&fnorm);
  fnorm_ = fnorm;

  // initial residual is r = f - M x for the right preconditioner
  // and r = H (f - M x)  for the left preconditioner
  if (left_pc_) {
    m_->Apply(x, p);
    p.Update(1.0, f, -1.0);
    h_->ApplyInverse(p, r);
  } else {
    m_->Apply(x, r);
    r.Update(1.0, f, -1.0);
  }

  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;
  if (rnorm0_ < 0.0) rnorm0_ = rnorm0;

  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Converged, itr=" << num_itrs_ << " ||r||=" << rnorm0 << std::endl;
    return criteria;  // Zero solution satifies all criteria.
  }

  // Ignore all criteria if one iteration is enforced.
  if (! (criteria & LIN_SOLVER_MAKE_ONE_ITERATION)) {
    int ierr = CheckConvergence_(rnorm0, fnorm);
    if (ierr != 0) return ierr;
  }

  // calculate the first (num_ritz_ + l) rows of c.
  int i0(0);
  double tmp;
  d.PutScalar(0.0);
  if (num_ritz_ == 0) {
    v_[0] = Teuchos::rcp(new Vector(r));
    v_[0]->Update(0.0, r, 1.0 / rnorm0);
    d(0) = rnorm0;
  } else {
    for (int i = 0 ; i <= num_ritz_; ++i) {
      r.Dot(*(v_[i]), &tmp);
      d(i) = tmp;
    }
    i0 = num_ritz_;
  }

  // set the leading diagonal block of T
  T.PutScalar(0.0);
  for (int i = 0; i <= num_ritz_; ++i) {
    for (int j = 0; j < num_ritz_; ++j) {
      T(i, j) = Hu_(i, j);
    }
  }

  // Apply Arnoldi method to extend the Krylov space calculate
  // at the end of the previous loop.
  double beta;
  for (int i = i0; i < krylov_dim_; i++) {
    if (left_pc_) {
      m_->Apply(*(v_[i]), p);
      h_->ApplyInverse(p, w);
    } else {
      h_->ApplyInverse(*(v_[i]), p);
      m_->Apply(p, w);
    }

    for (int k = 0; k <= i; k++) {  // Arnoldi algorithm
      w.Dot(*(v_[k]), &tmp);
      w.Update(-tmp, *(v_[k]), 1.0);
      T(k, i) = tmp;
    }
    w.Norm2(&beta);
    T(i + 1, i) = beta;
    if (beta == 0.0) break;

    v_[i + 1] = Teuchos::rcp(new Vector(w));
    v_[i + 1]->Update(0.0, r, 1.0 / beta);
  }

  // Solve the least-square problem min_d ||T d - c||.
  WhetStone::DenseMatrix Ttmp(T);
  int m(krylov_dim_ + 1), n(krylov_dim_), nrhs(1), info;
  int lwork(m * n);
  WhetStone::DenseVector work(lwork);

  WhetStone::DGELS_F77("N", &m, &n, &nrhs, Ttmp.Values(), &m, d.Values(), &m,
                       work.Values(), &lwork, &info);

  residual_ = fabs(d(n));
  num_itrs_ += krylov_dim_;
  ComputeSolution_(x, d.Values(), p, r); 

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << num_itrs_ << " ||r||=" << residual_ 
               << "  ritz vectors=" << num_ritz_ << std::endl;
  }
  int ierr = CheckConvergence_(residual_, fnorm);
  if (ierr != 0) return ierr;

  // Compute Schur vectors
  // -- allocate memory: Tm, Hm, and Vm
  WhetStone::DenseMatrix Tm(T, 0, krylov_dim_, 0, krylov_dim_);
  WhetStone::DenseMatrix Sm(Tm);
  WhetStone::DenseMatrix Vr(krylov_dim_ + 1, krylov_dim_);

  // -- auxiliary vector g = Tm^{-T} e_m
  Tm.Inverse();
  for (int i = 0; i < krylov_dim_; ++i) g(i) = beta * Tm(krylov_dim_ - 1, i);

  // -- solve eigenvector problem
  for (int i = 0; i < krylov_dim_; ++i) Sm(i, krylov_dim_ - 1) += beta * g(i);
  
  double Vl[1];
  WhetStone::DenseVector wr(krylov_dim_), wi(krylov_dim_);
  WhetStone::DGEEV_F77("N", "V", &n, Sm.Values(), &n, 
                       wr.Values(), wi.Values(), Vl, &nrhs, Vr.Values(), &m,
                       work.Values(), &lwork, &info);

  // -- select not more than (deflation_) Schur vectors and
  //    make them the first columns in Vr
  num_ritz_ = deflation_;

  for (int i = 0; i < num_ritz_; ++i) {
    int imin = i;
    double emin = wr(i);
    for (int j = i + 1; j < krylov_dim_; ++j) {
      if (wr(j) < emin) { 
        emin = wr(j);
        imin = j;
      }
    }
    wr.SwapRows(imin, i);
    wi.SwapRows(imin, i);

    Vr.SwapColumns(imin, i);
  }
  if (wi(num_ritz_ - 1) > 0.0) num_ritz_--;

  // -- add one vector and orthonormalize all columns.
  for (int i = 0; i < krylov_dim_; ++i) {
    Vr(krylov_dim_, i) = 0.0;
    Vr(i, num_ritz_) = -g(i);
  }
  Vr(krylov_dim_, num_ritz_) = 1.0;

  Vr.OrthonormalizeColumns(0, num_ritz_ + 1);

  // Calculate new basis for the next loop.
  std::vector<Teuchos::RCP<Vector> > vv(num_ritz_ + 1);
  for (int i = 0; i <= num_ritz_; ++i) {
    vv[i] = Teuchos::rcp(new Vector(x.Map()));
    vv[i]->PutScalar(0.0);
    for (int k = 0; k <= krylov_dim_; ++k) {
      vv[i]->Update(Vr(k, i), *(v_[k]), 1.0);
    }
  }
  
  for (int i = 0; i <= num_ritz_; ++i) {
    *(v_[i]) = *(vv[i]);
  }

  // Calculate modified Hessenberg matrix Hu = Vr_{nr+1}^T * T * Vr_nr
  WhetStone::DenseMatrix TVr(krylov_dim_ + 1, num_ritz_);
  WhetStone::DenseMatrix VTVr(num_ritz_ + 1, num_ritz_);

  WhetStone::DenseMatrix Vr1(Vr, 0, krylov_dim_, 0, num_ritz_);
  WhetStone::DenseMatrix Vr2(krylov_dim_ + 1, num_ritz_ + 1, Vr.Values(), WhetStone::WHETSTONE_DATA_ACCESS_VIEW);

  TVr.Multiply(T, Vr1, false);
  VTVr.Multiply(Vr2, TVr, true);
  Hu_ = VTVr;
  
  return LIN_SOLVER_MAX_ITERATIONS;
}


/* ******************************************************************
* Initialization from a parameter list. Available parameters:
* "error tolerance" [double] default = 1e-6
* "maximum number of iterations" [int] default = 100
* "convergence criteria" Array(string) default = "{relative rhs}"
****************************************************************** */
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
void IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::set_inverse_parameters(Teuchos::ParameterList& plist)
{
  InvIt::set_inverse_parameters(plist);

  controller_start_ = plist.get<int>("controller training start", 0);
  controller_end_ = plist.get<int>("controller training end", 3);
  controller_end_ = std::max(controller_end_, controller_start_ + 1);

  left_pc_ = (plist.get<std::string>("preconditioning strategy", "left") == "left");
  deflation_ = plist.get<int>("maximum size of deflation space", 0);
  num_ritz_ = 0;
}


/* ******************************************************************
* Givens rotations: initialization
****************************************************************** */
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
void IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::InitGivensRotation_(
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
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
void IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::ApplyGivensRotation_(
    double& dx, double& dy, double& cs, double& sn) const
{
  double tmp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = tmp;
}


/* ******************************************************************
* Computation of the solution destroys vector s.
* Right preconditioner uses two auxiliary vectors p and r.
****************************************************************** */
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
void IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::ComputeSolution_(
    Vector& x, int k, WhetStone::DenseMatrix& T, double* s,
    Vector& p, Vector& r) const
{
  for (int i = k; i >= 0; i--) {
    s[i] /= T(i, i);

    for (int j = i - 1; j >= 0; j--) {
      s[j] -= T(j, i) * s[i];
    }
  }

  // solution is x = x0 + V s for the left preconditioner
  //         and x = x0 + H V s for the right preconditioner
  if (left_pc_) {
    for (int j = 0; j <= k; j++) {
      x.Update(s[j], *(v_[j]), 1.0);
    }
  } else {
    p.PutScalar(0.0);
    for (int j = 0; j <= k; j++) {
      p.Update(s[j], *(v_[j]), 1.0);
    }
    h_->ApplyInverse(p, r);
    x.Update(1.0, r, 1.0);
  }
}


/* ******************************************************************
* solution is x = x0 + V s for the left preconditioner
*         and x = x0 + H V s for the right preconditioner
****************************************************************** */
template<class Matrix,class Preconditioner,class Vector,class VectorSpace>
void IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>::ComputeSolution_(
    Vector& x, double* d, Vector& p, Vector& r) const
{
  if (left_pc_) {
    for (int j = 0; j < krylov_dim_; j++) {
      x.Update(d[j], *(v_[j]), 1.0);
    }
  } else {
    p.PutScalar(0.0);
    for (int j = 0; j < krylov_dim_; j++) {
      p.Update(d[j], *(v_[j]), 1.0);
    }
    h_->ApplyInverse(p, r);
    x.Update(1.0, r, 1.0);
  }
}


}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
