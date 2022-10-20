#ifndef AMANZI_TEST_SOLVER_FNBASE1_HH_
#define AMANZI_TEST_SOLVER_FNBASE1_HH_

#include <math.h>
#include "Epetra_Vector.h"

#include "SolverFnBase.hh"

// ODE: f(u) = u (u^2 + 1) = 0.
class NonlinearProblem : public Amanzi::AmanziSolvers::SolverFnBase<Epetra_Vector> {
 public:
  NonlinearProblem(double atol, double rtol, bool exact_jacobian) :
    atol_(atol), rtol_(rtol), exact_jacobian_(exact_jacobian) {}

  void Residual(const Teuchos::RCP<Epetra_Vector>& u,
                const Teuchos::RCP<Epetra_Vector>& f) {
    for (int c = 0; c != u->MyLength(); ++c) {
      double x = (*u)[c];
      (*f)[c] = x * (x * x + 1.0);
    }
  }

  int ApplyPreconditioner(const Teuchos::RCP<const Epetra_Vector>& u,
                           const Teuchos::RCP<Epetra_Vector>& hu) {
    return hu->ReciprocalMultiply(1.0, *h_, *u, 0.0);
  }

  double ErrorNorm(const Amanzi::AmanziSolvers::Data<Epetra_Vector>& data,
                   const Amanzi::AmanziSolvers::ConvergenceMonitor& monitor) {
    auto u = data.u;
    auto r = data.r;

    double norm_r, norm_u;
    r->NormInf(&norm_r);
    u->NormInf(&norm_u);
    return norm_r;
    // return norm_r / (atol_ + rtol_ * norm_u);
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Epetra_Vector>& up) {
    h_ = Teuchos::rcp(new Epetra_Vector(*up));

    if (exact_jacobian_) {
      for (int c = 0; c != up->MyLength(); ++c) {
        double x = (*up)[c];
        (*h_)[c] = 3 * x * x + 1.0;
      }
    } else {
      for (int c = 0; c != up->MyLength(); ++c) {
        double x = (*up)[c];
        (*h_)[c] = x * x + 2.5;
      }
    }
  }

  void ChangedSolution() {};

 protected:
  double atol_, rtol_;
  bool exact_jacobian_;
  Teuchos::RCP<Epetra_Vector> h_;  // preconditioner
};

#endif

