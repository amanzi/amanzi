/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#ifndef AMANZI_TEST_SOLVER_FNBASE1_HH_
#define AMANZI_TEST_SOLVER_FNBASE1_HH_

#include <math.h>
#include "AmanziVector.hh"

#include "SolverFnBase.hh"

using namespace Amanzi;

// ODE: f(u) = u (u^2 + 1) = 0.
class NonlinearProblem : public AmanziSolvers::SolverFnBase<Vector_type> {
 public:
  NonlinearProblem(double atol, double rtol, bool exact_jacobian)
    : rtol_(rtol), atol_(atol), exact_jacobian_(exact_jacobian)
  {}

  void Residual(const Teuchos::RCP<Vector_type>& u,
                const Teuchos::RCP<Vector_type>& f)
  {
    auto uv = u->getLocalViewDevice();
    auto fv = f->getLocalViewDevice();
    Kokkos::parallel_for(
      "solver_fnbase1::Residual",
      fv.extent(0), KOKKOS_LAMBDA(const int c) {
      double x = uv(c, 0);
      fv(c, 0) = x * (x * x + 1.0);
    });
  }

  int ApplyPreconditioner(const Teuchos::RCP<const Vector_type>& u,
                          const Teuchos::RCP<Vector_type>& hu)
  {
    hu->elementWiseMultiply(1., *h_, *u, 0.);
    return 0;
  }

  double ErrorNorm(const Teuchos::RCP<const Vector_type>& u,
                   const Teuchos::RCP<const Vector_type>& du)
  {
    return du->normInf();
    // return norm_du / (atol_ + rtol_ * norm_u);
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Vector_type>& up)
  {
    if (!h_.get()) h_ = Teuchos::rcp(new Vector_type(up->getMap()));

    if (exact_jacobian_) {
      auto upv = up->getLocalViewDevice();
      auto hv = h_->getLocalViewDevice();
      Kokkos::parallel_for(
        "solver_fnbase1::UpdatePreconditioner loop 1",
        hv.extent(0), KOKKOS_LAMBDA(const int c) {
        double x = upv(c, 0);
        hv(c, 0) = 3 * x * x + 1.0;
      });
    } else {
      auto upv = up->getLocalViewDevice();
      auto hv = h_->getLocalViewDevice();
      Kokkos::parallel_for(
        "solver_fnbase1::UpdatePreconditioner loop 2",
        hv.extent(0), KOKKOS_LAMBDA(const int c) {
        double x = upv(c, 0);
        hv(c, 0) = x * x + 2.5;
      });
    }
    h_->reciprocal(*h_);
  }

  void ChangedSolution(){};

 protected:
  double atol_, rtol_;
  bool exact_jacobian_;
  Teuchos::RCP<Vector_type> h_; // preconditioner
};

#endif
