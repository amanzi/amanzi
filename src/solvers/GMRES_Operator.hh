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
#include "VerboseObject.hh"

#include "DenseMatrix.hh"
#include "LinearOperator.hh"
#include "LinearOperatorDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class GMRES_Operator : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  GMRES_Operator(Teuchos::RCP<const Matrix> m) :
      LinearOperator<Matrix, Vector, VectorSpace>(m) { 
    tol_ = 1e-6; 
    max_itrs_ = 100;
    krylov_dim_ = 10;
    criteria_ = SOLVER_CONVERGENCE_RELATIVE_RESIDUAL;
    initialized_ = false;
  }
  ~GMRES_Operator() {};

  void Init(Teuchos::ParameterList& plist);  

  void ApplyInverse(const Vector& v, Vector& hv) const { 
    num_itrs_ = gmres_restart(v, hv, tol_, max_itrs_, criteria_); 
  }

  Teuchos::RCP<GMRES_Operator> Clone() const {};

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void set_krylov_dim(int n) { krylov_dim_ = n; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }

 public:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  int gmres_restart(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  int gmres(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;
  void ComputeSolution(Vector& x, int k, WhetStone::DenseMatrix& T, double* s, Vector** v) const;
  void InitGivensRotation( double& dx, double& dy, double& cs, double& sn) const;
  void ApplyGivensRotation(double& dx, double& dy, double& cs, double& sn) const;

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  int max_itrs_, criteria_, krylov_dim_;
  double tol_;
  mutable int num_itrs_;
  mutable double residual_;
  mutable bool initialized_;
};


/* ******************************************************************
* GMRES with restart input/output data:
*  f [input]         the right-hand side
*  x [input/output]  initial guess / final solution
*  tol [input]       convergence tolerance
*  max_itrs [input]  maximum number of iterations
*  criteria [input]  sum of termination critaria
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
int GMRES_Operator<Matrix, Vector, VectorSpace>::gmres_restart(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  int num_itrs = 0;
  int max_itrs_left = max_itrs;

  int i = krylov_dim_;
  while (i == krylov_dim_ && max_itrs_left > 0) {
    i = gmres(f, x, tol, max_itrs_left, criteria); 
    num_itrs += i;
    max_itrs_left -= i;
  }

  return num_itrs;
}


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
  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;
  
  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    return 0;
  }
  if (! (criteria & SOLVER_MAKE_ONE_ITERATION)) {
    if (criteria & SOLVER_CONVERGENCE_RELATIVE_RHS) {
      if (rnorm0 < tol * fnorm) return 0; 
    } else if (criteria & SOLVER_CONVERGENCE_ABSOLUTE_RESIDUAL) {
      if (rnorm0 < tol) return 0;
    }
  }

  v[0] = new Vector(r);
  v[0]->Update(0.0, r, 1.0 / rnorm0);

  s[0] = rnorm0;
    
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
    s[i + 1] = 0.0;

    v[i + 1] = new Vector(w);
    v[i + 1]->Update(0.0, r, 1.0 / tmp);

    for (int k = 0; k < i; k++) {
      ApplyGivensRotation(T(k, i), T(k + 1, i), cs[k], sn[k]);
    }
      
    InitGivensRotation( T(i, i), T(i + 1, i), cs[i], sn[i]);
    ApplyGivensRotation(T(i, i), T(i + 1, i), cs[i], sn[i]);
    ApplyGivensRotation(s[i],    s[i + 1],    cs[i], sn[i]);
    residual_ = fabs(s[i + 1]);

    if (initialized_) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << i << " ||r||=" << residual_ << endl;
      }
    }
    if (criteria_ & SOLVER_CONVERGENCE_RELATIVE_RHS) {
      if (residual_ < tol * fnorm) { 
        ComputeSolution(x, i, T, s, v);  // vector s is overwritten
        return i;
      }
    } else if (criteria_ & SOLVER_CONVERGENCE_RELATIVE_RESIDUAL) {
      if (residual_ < tol * rnorm0) { 
        ComputeSolution(x, i, T, s, v);  // vector s is overwritten
        return i;
      }
    } else if (criteria_ & SOLVER_CONVERGENCE_ABSOLUTE_RESIDUAL) {
      if (residual_ < tol) { 
        ComputeSolution(x, i, T, s, v);  // vector s is overwritten
        return i;
      }
    }
  }

  ComputeSolution(x, krylov_dim_ - 1, T, s, v);  // vector s is overwritten
  return krylov_dim_;
}


/* ******************************************************************
* Initialization from a parameter list. Available parameters:
* "error tolerance" [double] default = 1e-6
* "maximum number of iterations" [int] default = 100
* "convergence criteria" Array(string) default = "{relative rhs}"
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void GMRES_Operator<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Amanzi::PCG_Solver", plist)); 

  double tol = plist.get<double>("error tolerance", 1e-6);
  set_tolerance(tol);

  double max_itrs = plist.get<int>("maximum number of iterations", 100);
  set_max_itrs(max_itrs);

  double krylov_dim = plist.get<int>("size of Krylov space", 10);
  set_krylov_dim(krylov_dim);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names = plist.get<Teuchos::Array<std::string> > ("convergence criteria").toVector();

    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "relative rhs") {
        criteria += SOLVER_CONVERGENCE_RELATIVE_RHS;
      } else if (names[i] == "relative residual") {
        criteria += SOLVER_CONVERGENCE_RELATIVE_RESIDUAL;
      } else if (names[i] == "absolute residual") {
        criteria += SOLVER_CONVERGENCE_ABSOLUTE_RESIDUAL;
      } else if (names[i] == "make one iteration") {
        criteria += SOLVER_MAKE_ONE_ITERATION;
      }
    }
  } else {
    criteria = SOLVER_CONVERGENCE_RELATIVE_RESIDUAL;
  }

  set_criteria(criteria);
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
* Computation of the solution destroys vector s.
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void GMRES_Operator<Matrix, Vector, VectorSpace>::ComputeSolution(
    Vector& x, int k, WhetStone::DenseMatrix& T, double* s, Vector** v) const
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

