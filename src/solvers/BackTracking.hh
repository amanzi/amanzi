/*
  This is the Nonlinear Solver component of the Amanzi code.

  Interface for backtracking algorithms.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#ifndef AMANZI_BACKTRACKING_
#define AMANZI_BACKTRACKING_

#include "Teuchos_RCP.hpp"

#include "SolverDefs.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector>
class BackTracking {
 public:
  BackTracking(const Teuchos::RCP<SolverFnBase<Vector> >& fn)
      : fn_(fn) {};

  int Bisection(const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du);
  int Bisection(double f0, const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du);

  // access
  int num_steps() { return num_steps_; }
  double initial_residual() { return initial_residual_; }
  double final_residual() { return final_residual_; }

 private:
  Teuchos::RCP<SolverFnBase<Vector> > fn_;

  int num_steps_;
  double initial_residual_, final_residual_;
};


/* ******************************************************************
* Bisection algorithm: u1 = u0 - s * du
****************************************************************** */
template<class Vector>
int BackTracking<Vector>::Bisection(double f0, const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du)
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
    return 0;  // success
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
      return BACKTRACKING_USED;  // backtraking was activated
    } else {
      s1 = s2;
    }
  }  

  return BACKTRACKING_MAX_ITERATIONS;
}


/* ******************************************************************
* Bisection algorithm
****************************************************************** */
template<class Vector>
int BackTracking<Vector>::Bisection(const Teuchos::RCP<Vector> u0, Teuchos::RCP<Vector> du)
{
  double f0;
  Teuchos::RCP<Vector> r = Teuchos::rcp(new Vector(*u0));

  fn_->Residual(u0, r);
  r->Norm2(&f0);

  return Bisection(f0, u0, du);  
}
  
}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
