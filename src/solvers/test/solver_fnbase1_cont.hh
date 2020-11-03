#ifndef AMANZI_TEST_SOLVER_FNBASE1_HH_
#define AMANZI_TEST_SOLVER_FNBASE1_HH_

#include <math.h>
#include "AmanziVector.hh"

#include "SolverFnBase.hh"

using namespace Amanzi; 

// ODE: f(u) = u (u^2 + 1) = 0.
class NonlinearProblem : public Amanzi::AmanziSolvers::SolverFnBase<Vector_type> {
 public:
  NonlinearProblem(double atol, double rtol, bool exact_jacobian) :
    rtol_(rtol), atol_(atol), exact_jacobian_(exact_jacobian) {}

  void Residual(const Teuchos::RCP<Vector_type>& u,
                const Teuchos::RCP<Vector_type>& f) {

    auto uv = u->getLocalViewHost();
    auto fv = f->getLocalViewHost(); 

    for (int c = 0; c != u->getLocalLength(); ++c) {
      double x = uv(c,0);
      fv(c,0) = std::pow(x,lambda_) - 0.5;
      std::cout << "x=" << x << ", f=" << fv(c,0) << std::endl;
    }
  }

  int ApplyPreconditioner(const Teuchos::RCP<const Vector_type>& u,
                           const Teuchos::RCP<Vector_type>& hu) {
    return hu->ReciprocalMultiply(1.0, *h_, *u, 0.0);
  }

  double ErrorNorm(const Teuchos::RCP<const Vector_type>& u,
                   const Teuchos::RCP<const Vector_type>& du) {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    u->NormInf(&norm_u);
    return norm_du / (atol_ + rtol_ * norm_u);
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Vector_type>& up) {
    h_ = Teuchos::rcp(new Vector_type(*up));

    if (exact_jacobian_) {
      for (int c = 0; c != up->getLocalLength(); ++c) {
        double x = (*up)[c];
        (*h_)[c] = lambda_ * std::pow(x,lambda_-1);
      }
    } else {
      for (int c = 0; c != up->getLocalLength(); ++c) {
        (*h_)[c] = 1.0;
      }
    }
  }

  void UpdateContinuationParameter(double lambda) {
    lambda_ = 3. - 2.*lambda;
    std::cout << "Updating lambda = " << lambda_ << std::endl;
  }

  void ChangedSolution() {};

 protected:
  double atol_, rtol_;
  bool exact_jacobian_;
  double lambda_;
  Teuchos::RCP<Vector_type> h_;  // preconditioner
};

#endif

