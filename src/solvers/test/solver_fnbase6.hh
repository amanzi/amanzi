#ifndef AMANZI_TEST_SOLVER_FNBASE6_HH_
#define AMANZI_TEST_SOLVER_FNBASE6_HH_

#include <math.h>
#include "Epetra_Vector.h"

#include "SolverFnBase.hh"

// ODE: f(u) = ... = 0
class NonlinearProblem6 : public Amanzi::AmanziSolvers::SolverFnBase<Epetra_Vector> {
 public:
  NonlinearProblem6(double atol, double rtol, bool exact_jacobian)
    : atol_(atol), rtol_(rtol), exact_jacobian_(exact_jacobian)
  {}

  void Residual(const Teuchos::RCP<Epetra_Vector>& u, const Teuchos::RCP<Epetra_Vector>& f)
  {
    for (int c = 0; c != u->MyLength(); ++c) {
      double x = (*u)[c];
      (*f)[c] = x < 0 ? -pow(fabs(x), 0.2) : pow(fabs(x), 0.2);
    }
    std::cout << "  Evaluating: f(u=" << (*u)[0] << "," << (*u)[1] << ") = " << (*f)[0] << ","
              << (*f)[1] << std::endl;
  }

  int ApplyPreconditioner(const Teuchos::RCP<const Epetra_Vector>& u,
                          const Teuchos::RCP<Epetra_Vector>& hu)
  {
    hu->ReciprocalMultiply(1.0, *h_, *u, 0.0);
    return 0;
  }

  double
  ErrorNorm(const Teuchos::RCP<const Epetra_Vector>& u, const Teuchos::RCP<const Epetra_Vector>& du)
  {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    double error = norm_du / atol_;
    std::cout << " : error = " << error << std::endl;
    return error;
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Epetra_Vector>& up)
  {
    h_ = Teuchos::rcp(new Epetra_Vector(*up));

    if (exact_jacobian_) {
      for (int c = 0; c != up->MyLength(); ++c) {
        double x = (*up)[c];
        (*h_)[c] = 0.2 / pow(fabs(x), 0.8);
      }
    } else {
      for (int c = 0; c != up->MyLength(); ++c) {
        double x = (*up)[c];
        (*h_)[c] = 0.3 / pow(fabs(x), 2.0 / 3);
      }
    }
  }

  void ChangedSolution(){};

 protected:
  double atol_, rtol_;
  bool exact_jacobian_;
  Teuchos::RCP<Epetra_Vector> h_; // preconditioner
};

#endif
