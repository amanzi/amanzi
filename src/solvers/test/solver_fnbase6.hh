#ifndef AMANZI_TEST_SOLVER_FNBASE6_HH_
#define AMANZI_TEST_SOLVER_FNBASE6_HH_

#include <math.h>
#include "AmanziVector.hh"

#include "SolverFnBase.hh"

// ODE: f(u) = ... = 0
class NonlinearProblem6 : public Amanzi::AmanziSolvers::SolverFnBase<Vector_type> {
 public:
  NonlinearProblem6(double atol, double rtol, bool exact_jacobian) :
    rtol_(rtol), atol_(atol), exact_jacobian_(exact_jacobian) {}

  void Residual(const Teuchos::RCP<Vector_type>& u,
                const Teuchos::RCP<Vector_type>& f) {
    auto uv = u->getLocalViewDevice();
    auto fv = f->getLocalViewDevice();
    Kokkos::parallel_for(
      "solver_fnbase6::Residual",
      fv.extent(0), KOKKOS_LAMBDA(const int c) {
      double x = uv(c, 0);
      fv(c, 0) = x < 0 ? -pow(fabs(x), 0.2) : pow(fabs(x), 0.2);
    });
  }

  int ApplyPreconditioner(const Teuchos::RCP<const Vector_type>& u,
                           const Teuchos::RCP<Vector_type>& hu) {
    hu->elementWiseMultiply(1.0, *h_, *u, 0.0);
    return 0;
  }

  double ErrorNorm(const Teuchos::RCP<const Vector_type>& u,
                   const Teuchos::RCP<const Vector_type>& du) {
    double norm_du;
    norm_du = du->normInf();
    double error = norm_du / atol_;
    std::cout << " : error = " << error << std::endl;
    return error;
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Vector_type>& up) {
    if (!h_.get()) h_ = Teuchos::rcp(new Vector_type(up->getMap()));

    if (exact_jacobian_) {
      auto upv = up->getLocalViewDevice();
      auto hv = h_->getLocalViewDevice();
      Kokkos::parallel_for(
      "solver_fnbase6::UpdatePreconditioner loop 1",
        hv.extent(0), KOKKOS_LAMBDA(const int c) {
        double x = upv(c, 0);
        hv(c, 0) = alpha_ * pow(fabs(x), alpha_ - 1.);
      });
    } else {
      auto upv = up->getLocalViewDevice();
      auto hv = h_->getLocalViewDevice();
      Kokkos::parallel_for(
      "solver_fnbase6::UpdatePreconditioner loop 2",
        hv.extent(0), KOKKOS_LAMBDA(const int c) {
        double x = upv(c, 0);
        hv(c, 0) = 0.3 * pow(fabs(x), -.6666667);
      });
    }
    h_->reciprocal(*h_);
  }

  void ChangedSolution() {};

 protected:
  double atol_, rtol_;
  bool exact_jacobian_;
  Teuchos::RCP<Vector_type> h_;  // preconditioner
  
  static constexpr double alpha_ = 0.2; 

};

#endif

