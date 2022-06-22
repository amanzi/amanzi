#ifndef AMANZI_TEST_SOLVER_FNBASE1_HH_
#define AMANZI_TEST_SOLVER_FNBASE1_HH_

#include <math.h>
#include "AmanziVector.hh"
#include "AmanziDebug.hh"
#include "SolverFnBase.hh"

using namespace Amanzi; 

// ODE: f(u) = u (u^2 + 1) = 0.
class NonlinearProblem : public Amanzi::AmanziSolvers::SolverFnBase<Vector_type> {
 public:
  NonlinearProblem(double atol, double rtol, bool exact_jacobian) :
    rtol_(rtol), atol_(atol), exact_jacobian_(exact_jacobian) {}

  void Residual(const Teuchos::RCP<Vector_type>& u,
                const Teuchos::RCP<Vector_type>& f) {
    auto uv = u->getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto fv = f->getLocalViewDevice(Tpetra::Access::ReadWrite);
    assert(uv != fv);
    Kokkos::parallel_for(
      "solver_fnbase1::Residual",
      fv.extent(0), KOKKOS_LAMBDA(const int c) {
      double x = uv(c, 0);
      fv(c, 0) = x * (x * x + 1.0);
    });
  }

  int ApplyPreconditioner(const Teuchos::RCP<const Vector_type>& u,
                           const Teuchos::RCP<Vector_type>& hu) {
    std::cout << std::setprecision(16) << "SolverFnBase1: u = " << Debug::get0(*u) << std::endl;
    hu->elementWiseMultiply(1., *h_, *u, 0.);
    std::cout << "SolverFnBase1: hu = " << Debug::get0(*hu) << std::endl;
    return 0;
  }

  double ErrorNorm(const Teuchos::RCP<const Vector_type>& u,
                   const Teuchos::RCP<const Vector_type>& du) {
    return du->normInf();
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Vector_type>& up) {
    if (!h_.get()) h_ = Teuchos::rcp(new Vector_type(up->getMap()));
    std::cout << "SolverFnBase1: up = " << Amanzi::Debug::get0(*up) << std::endl;

    if (exact_jacobian_) {
      auto upv = up->getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto hv = h_->getLocalViewDevice(Tpetra::Access::ReadWrite);
      Kokkos::parallel_for(
        "solver_fnbase1::UpdatePreconditioner loop 1",
        hv.extent(0), KOKKOS_LAMBDA(const int c) {
        double x = upv(c, 0);
        hv(c, 0) = 3 * x * x + 1.0;
      });
    } else {
      auto upv = up->getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto hv = h_->getLocalViewDevice(Tpetra::Access::ReadWrite);
      Kokkos::parallel_for(
        "solver_fnbase1::UpdatePreconditioner loop 2",
        hv.extent(0), KOKKOS_LAMBDA(const int c) {
        double x = upv(c, 0);
        hv(c, 0) = x * x + 2.5;
      });
    }
    std::cout << "SolverFnBase1: h = " << Amanzi::Debug::get0(*h_) << std::endl;
    h_->reciprocal(*h_);
  }

  void ChangedSolution() {};

 protected:
  double atol_, rtol_;
  bool exact_jacobian_;
  Teuchos::RCP<Vector_type> h_;  // preconditioner
};

#endif

