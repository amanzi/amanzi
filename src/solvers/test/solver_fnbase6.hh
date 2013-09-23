#ifndef AMANZI_TEST_SOLVER_FNBASE2_HH_
#define AMANZI_TEST_SOLVER_FNBASE2_HH_

#include <math.h>
#include "Epetra_Vector.h"

#include "SolverFnBase.hh"

// ODE for testing
class NonlinearProblem : public Amanzi::AmanziSolvers::SolverFnBase<Epetra_Vector> {
 public:
  NonlinearProblem(double atol, double rtol, bool real_precon) :
    rtol_(rtol), atol_(atol), real_precon_(real_precon) {}

  void Residual(const Teuchos::RCP<Epetra_Vector>& u,
                const Teuchos::RCP<Epetra_Vector>& f) {
    // f = tanh(u)
    for (int c = 0; c != u->MyLength(); ++c) {
      double x = (*u)[c];
      (*f)[c] = x < 0 ? -std::pow(std::abs(x), 1./5) : std::pow(std::abs(x), 1./5);
    }

    // std::cout << " Residual eval:" << std::endl;
    // std::cout << "  f(" << (*u)[0] << ") = " << (*f)[0] << std::endl;
    // std::cout << "  f(" << (*u)[1] << ") = " << (*f)[1] << std::endl;
  }

  void ApplyPreconditioner(const Teuchos::RCP<const Epetra_Vector>& u,
                           const Teuchos::RCP<Epetra_Vector>& Pu) {
    Pu->ReciprocalMultiply(1.0, *Pu_, *u, 0.0);
  }

  double ErrorNorm(const Teuchos::RCP<const Epetra_Vector>& u,
                   const Teuchos::RCP<const Epetra_Vector>& du) {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    return norm_du / (atol_ + rtol_ * norm_u);
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Epetra_Vector>& up) {
    if (Pu_ == Teuchos::null) {
      Pu_ = Teuchos::rcp(new Epetra_Vector(*up));
    }

    if (real_precon_) {
      for (int c = 0; c != up->MyLength(); ++c) {
        double x = (*up)[c];
        (*Pu_)[c] = 1.0 / std::pow(std::abs(x), 4.0 / 5) / 5.0;
        // std::cout << "Pu = " << (*Pu_)[c] << std::endl;
      }
    } else {
      for (int c = 0; c != up->MyLength(); ++c) {
        double x = (*up)[c];
        (*Pu_)[c] = 1.0 / std::pow(std::abs(x), 2.0 / 3) / 3.0;
      }
    }
  }

  bool ModifyCorrection(const Teuchos::RCP<const Epetra_Vector>& res,
                        const Teuchos::RCP<const Epetra_Vector>& u,
			const Teuchos::RCP<Epetra_Vector>& du) {
    return false; 
  }

  void ChangedSolution() {};

  double atol_, rtol_;
  bool real_precon_;  // use "approximate" Jacobian
  Teuchos::RCP<Epetra_Vector> Pu_;
};

#endif

