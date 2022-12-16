/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Solvers

  Interface for backtracking algorithms.
*/


#ifndef AMANZI_BACKTRACKING_
#define AMANZI_BACKTRACKING_

#include "Teuchos_RCP.hpp"

#include "SolverDefs.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector>
class BackTracking {
 public:
  BackTracking(const Teuchos::RCP<SolverFnBase<Vector>>& fn) : fn_(fn){};

  int Bisection(const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du);
  int Bisection(double f0, const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du);

  int LineSearch(Vector& xold,
                 double fold,
                 Vector& g,
                 Vector& p,
                 Vector& x,
                 double& f,
                 double step_max);

  // access
  int num_steps() { return num_steps_; }
  double initial_residual() { return initial_residual_; }
  double final_residual() { return final_residual_; }
  double fun_calls() { return fun_calls_; }

 private:
  Teuchos::RCP<SolverFnBase<Vector>> fn_;

  int num_steps_, fun_calls_;
  double initial_residual_, final_residual_;
};


/* ******************************************************************
* Bisection algorithm: u1 = u0 - s * du
****************************************************************** */
template <class Vector>
int
BackTracking<Vector>::Bisection(double f0, const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du)
{
  double f1;
  initial_residual_ = f0;
  final_residual_ = f0;

  Teuchos::RCP<Vector> u1 = Teuchos::rcp(new Vector(*u0));
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(*u0));

  u1->Update(1.0, *u0, -1.0, *du, 0.0);

  fn_->Residual(u1, r);
  r->Norm2(&f1);

  num_steps_ = 0;
  if (f1 < f0 * BACKTRACKING_GOOD_REDUCTION) {
    final_residual_ = f1;
    return 0; // success
  }

  double s0(0.0), s1(1.0), s2;
  for (int n = 0; n < BACKTRACKING_MAX_ITERATIONS; n++) {
    num_steps_++;
    s2 = (s0 + s1) / 2;
    u1->Update(1.0, *u0, -s2, *du, 0.0);

    fn_->Residual(u1, r);
    r->Norm2(&f1);

    if (f1 < f0 * BACKTRACKING_GOOD_REDUCTION) {
      du->Scale(s2);
      final_residual_ = f1;
      return BACKTRACKING_USED; // backtraking was activated
    } else {
      s1 = s2;
    }
  }

  return BACKTRACKING_MAX_ITERATIONS;
}


/* ******************************************************************
* Bisection algorithm
****************************************************************** */
template <class Vector>
int
BackTracking<Vector>::Bisection(const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du)
{
  double f0;
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(*u0));

  fn_->Residual(u0, r);
  r->Norm2(&f0);

  return Bisection(f0, u0, du);
}


/* ******************************************************************
* Line search
****************************************************************** */
template <class Vector>
int
BackTracking<Vector>::LineSearch(Vector& xold,
                                 double fold,
                                 Vector& g,
                                 Vector& p,
                                 Vector& x,
                                 double& f,
                                 double step_max)
{
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(x));

  // alpha ensures sufficient decrease in function value
  // tolx is the convergence criterion on x.
  double alpha = 1.0e-4, tolx = 1e-6;
  double tmplam, disc, alam, alam2 = 0.0, alamin, f2 = 0.0;

  double sum;
  p.Norm2(&sum);

  // Scale if attempted step is too big.
  if (sum > step_max) { p.Scale(step_max / sum); }

  double slope;
  g.Dot(p, &slope);
  if (slope >= 0.0) return BACKTRACKING_ROUNDOFF_PROBLEM;

  // Compute lambda_min
  /*
  double temp, test = 0.0;
  for (int i = 0; i < n; i++) {
    temp=abs(p[i])/MAX(abs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  */

  // Always try full Newton step first.
  // alamin = tolx / test;
  alamin = tolx;
  alam = 1.0;
  for (int n = 0; n < BACKTRACKING_MAX_ITERATIONS; n++) {
    x.Update(1.0, xold, alam, p, 0.0);
    fn_->Residual(x, r);
    r->Norm2(&f);

    // Convergence on Delta x.
    if (alam < alamin) {
      x = xold;
      return 0;
      // Sufficient function decrease
    } else if (f <= fold + alpha * alam * slope) {
      return 0;
      // Backtracking
    } else {
      if (alam == 1.0) {
        tmplam = -slope / (2.0 * (f - fold - slope)); // First time.
      } else {                                        // Subsequent backtracks.
        double rhs1, rhs2, a, b;
        rhs1 = f - fold - alam * slope;
        rhs2 = f2 - fold - alam2 * slope;
        a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
        b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);

        if (a == 0.0) {
          tmplam = -slope / (2.0 * b);
        } else {
          disc = b * b - 3.0 * a * slope;
          if (disc < 0.0) {
            tmplam = 0.5 * alam;
          } else if (b <= 0.0) {
            tmplam = (-b + sqrt(disc)) / (3.0 * a);
          } else {
            tmplam = -slope / (b + sqrt(disc));
          }
        }

        // lambda <= lambda_1 / 2
        if (tmplam > 0.5 * alam) tmplam = 0.5 * alam;
      }
    }

    alam2 = alam;
    f2 = f;
    alam = std::max(tmplam, 0.1 * alam); // lambda > 0.1 lambda_1
  }
  return 0;
}


} // namespace AmanziSolvers
} // namespace Amanzi

#endif
