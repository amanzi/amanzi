/*
This is the Linear Solver component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
         Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#ifndef __PCG_OPERATOR_HH__
#define __PCG_OPERATOR_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "errors.hh"
#include "Solver_constants.hh"
 
namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class PCG_Operator : public Matrix {
 public:
  PCG_Operator(Teuchos::RCP<const Matrix> m) : m_(m) { 
    tol_ = 1e-6; 
    max_itrs_ = 100;
    criteria_ = SOLVER_CONVERGENCE_RHS;
  }
  ~PCG_Operator() {};

  void Init(Teuchos::ParameterList& plist);  

  void Apply(const Vector& v, Vector& mv) const { m_->Apply(v, mv); }
  void ApplyInverse(const Vector& v, Vector& hv) const { 
    num_itrs_ = pcg(v, hv, tol_, max_itrs_, criteria_); 
  }

  Teuchos::RCP<const VectorSpace> domain() const { return m_->domain(); }
  Teuchos::RCP<const VectorSpace> range() const { return m_->range(); }
  Teuchos::RCP<PCG_Operator> Clone() const {};

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }

 private:
  int pcg(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

 private:
  Teuchos::RCP<const Matrix> m_;

  int max_itrs_, criteria_;
  double tol_;
  mutable int num_itrs_;
  mutable double residual_;
};


/* ******************************************************************
* PCG input/output data:
*  f [input]         the right-hand side
*  x [input/output]  initial guess / final solution
*  tol [input]       convergence tolerance
*  max_itrs [input]  maximum number of iterations
*  criteria [input]  sum of termination critaria
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
int PCG_Operator<Matrix, Vector, VectorSpace>::pcg(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Vector r(f), p(f), v(f);  // construct empty vectors

  double fnorm, xnorm;
  f.Norm2(&fnorm);
  x.Norm2(&xnorm);
  if (xnorm == 0) xnorm = fnorm;
  if (xnorm == 0) return 0;

  m_->Apply(x, r);  // r = f - M * x
  r.Update(1.0, f, -1.0);
  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;

  m_->ApplyInverse(r, p);  // gamma = (H r,r)
  double gamma0;
  p.Dot(r, &gamma0);
  if (gamma0 < 0) {
    Errors::Message msg("PCG: ApplyInverse() is not an SPD operator.");
    Exceptions::amanzi_throw(msg);
  }

  if (criteria & SOLVER_CONVERGENCE_RHS)   
    if( rnorm0 < tol * fnorm) return 0; 

  for (int i = 0; i < max_itrs; i++) {
    m_->Apply(p, v);
    double alpha;
    v.Dot(p, &alpha);

    if (alpha < 0.0) {
      Errors::Message msg("PCG: Apply() is not an SPD operator.");
      Exceptions::amanzi_throw(msg);
    }
    alpha = gamma0 / alpha;   

    x.Update( alpha, p, 1.0);
    r.Update(-alpha, v, 1.0);

    m_->ApplyInverse(r, v);  // gamma = (H r, r)
    double gamma1;
    v.Dot(r, &gamma1);
    if (gamma1 < 0.0) {
      Errors::Message msg("PCG: ApplyInverse() is not an SPD operator.");
      Exceptions::amanzi_throw(msg);
    }

    double rnorm;
    r.Norm2(&rnorm);
    residual_ = rnorm;
 
    if (criteria & SOLVER_CONVERGENCE_RHS) 
      if (rnorm < tol * fnorm) return i+1;
    if (criteria & SOLVER_CONVERGENCE_RESIDUAL) 
      if (rnorm < tol * rnorm0) return i+1;

    double beta = gamma1 / gamma0;
    gamma0 = gamma1;
 
    p.Update(1.0, v, beta);
  }

  return max_itrs;
};


/* ******************************************************************
* Initialization from a parameter list. Available parameters:
* "error tolerance" [double] default = 1e-6
* "maximum number of iterations" [int] default = 100
* "convergence criteria" Array(string) default = "{relative rhs}"
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void PCG_Operator<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  double tol = plist.get<double>("error tolerance", 1e-6);
  set_tolerance(tol);

  double max_itrs = plist.get<int>("maximum number of iterations", 100);
  set_max_itrs(max_itrs);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names = plist.get<Teuchos::Array<std::string> > ("convergence criteria").toVector();

    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "relative rhs") {
        criteria += SOLVER_CONVERGENCE_RHS;
      } else if (names[i] == "relative residual") {
        criteria += SOLVER_CONVERGENCE_RESIDUAL;
      }
    }
  } else {
    criteria = SOLVER_CONVERGENCE_RHS;
  }

  set_criteria(criteria);
}

}  // namespace AmanziSolvers
}  // namespace Amanzi
 
#endif
               

