/*
This is the Linear Solver component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
         Konstantin Lipnikov (lipnikov@lanl.gov)

Generalized minimum residual method.
Usage: 
*/

#ifndef  __GMRES_OPERATOR_HH__
#define  __GMRES_OPERATOR_HH__

#include <cmath>

#include "Teuchos_RCP.hpp"
#include "exceptions.hh"
#include "DenseMatrix.hh"

#include "Solver_constants.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class GMRES_Operator : public Matrix {
 public:
  GMRES_Operator(Teuchos::RCP<const Matrix> m) : m_(m) { 
    tol_ = 1e-6; 
    max_itrs_ = 100;
    krylov_dim_ = 10;
    criteria_ = SOLVER_CONVERGENCE_RHS;
  }
  ~GMRES_Operator() {};

  void Apply(const Vector& v, Vector& mv) const { m_->Apply(v, mv); }
  void ApplyInverse(const Vector& v, Vector& hv) const { 
    num_itrs_ = gmres(v, hv, tol_, max_itrs_, criteria_); 
  }

  Teuchos::RCP<const VectorSpace> domain() const { return m_->domain(); }
  Teuchos::RCP<const VectorSpace> range() const { return m_->range(); }
  Teuchos::RCP<GMRES_Operator> Clone() const {};

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void set_krylov_dim(int n) { krylov_dim_ = n; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }

 private:
  int gmres(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  void UpdateSolution(Vector& x, int k, WhetStone::DenseMatrix& T, double* s, Vector** v) const;
  void InitGivensRotation( double& dx, double& dy, double& cs, double& sn) const;
  void ApplyGivensRotation(double& dx, double& dy, double& cs, double& sn) const;

 private:
  Teuchos::RCP<const Matrix> m_;

  int max_itrs_, criteria_, krylov_dim_;
  double tol_;
  mutable int num_itrs_;
  mutable double residual_;
};


/* ******************************************************************
* GMRES input/output data:
*  f [input]         the right-hand side
*  x [input/output]  initial guess / final solution
*  tol [input]       convergence tolerance
*  max_itrs [input]  maximum number of iterations
*  criteria [input]  sum of termination critaria
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
int GMRES_Operator<Matrix, Vector, VectorSpace>::gmres(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Vector w(f), r(f), p(f);  // construct empty vectors
  Vector* v[krylov_dim_ + 1];

  double s[krylov_dim_ + 1], cs[krylov_dim_ + 1], sn[krylov_dim_ + 1];
  WhetStone::DenseMatrix T(krylov_dim_ + 1, krylov_dim_);

  m_->ApplyInverse(f, r);
  double fnorm;
  r.Norm2(&fnorm);
    
  m_->Apply(x, p);  // p = f - M * x
  p.Update(1.0, f, -1.0);

  m_->ApplyInverse(p, r);
  double beta;
  r.Norm2(&beta);
  
  if (fabs(fnorm) < 1.0e-30) fnorm = 1.0;
  if ((residual_ = beta / fnorm) <= tol_) return 0;

  v[0] = new Vector(r);
  v[0]->Update(0.0, r, 1.0 / beta);

  s[0] = beta;
    
  for (int i = 0; i < krylov_dim_; i++) {
    m_->Apply(*(v[i]), p);
    m_->ApplyInverse(p, w);

    double tmp; 
    for (int k = 0; k <= i; k++) {
      w.Dot(*(v[k]), &tmp);
      w.Update(-tmp, *(v[k]), 1.0);
      T(k, i) = tmp;
    }
    w.Norm2(&tmp);
    T(i + 1, i) = tmp;

    v[i + 1] = new Vector(w);
    v[i + 1]->Update(0.0, r, 1.0 / tmp);

    for (int k = 0; k < i; k++) {
      ApplyGivensRotation(T(k, i), T(k + 1, i), cs[k], sn[k]);
    }
      
    InitGivensRotation( T(i, i), T(i + 1, i), cs[i], sn[i]);
    ApplyGivensRotation(T(i, i), T(i + 1, i), cs[i], sn[i]);
    ApplyGivensRotation(s[i],    s[i + 1],    cs[i], sn[i]);

    // cout << "itr = " << i+1 << "  residual = " << fabs(s(i+1)) << endl;
    
    if (criteria_ & SOLVER_CONVERGENCE_RHS || 
        criteria_ & SOLVER_CONVERGENCE_RESIDUAL) {
      if ((residual_ = fabs(s[i + 1]) / fnorm) < tol_) { 
        UpdateSolution(x, i, T, s, v);
        return i;
      }
    }
  }

  UpdateSolution(x, krylov_dim_ - 1, T, s, v);
  return krylov_dim_;
}



/* ******************************************************************
* Givens rotations: initialization
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void GMRES_Operator<Matrix, Vector, VectorSpace>::InitGivensRotation(
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
void GMRES_Operator<Matrix, Vector, VectorSpace>::ApplyGivensRotation(
    double& dx, double& dy, double& cs, double& sn) const
{
  double tmp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = tmp;
}


/* ******************************************************************
* Backward solve
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void GMRES_Operator<Matrix, Vector, VectorSpace>::UpdateSolution(
    Vector& x, int k, WhetStone::DenseMatrix& T, double* s, Vector** v) const
{
  double y[k + 1];

  for (int i = k; i >= 0; i--) {
    y[i] = s[i] / T(i, i);

    for (int j = i - 1; j >= 0; j--) {
      y[j] -= T(j, i) * y[i];
    }
  }

  for (int j = 0; j <= k; j++) {
    x.Update(y[j], *(v[j]), 1.0);
  }
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif  

