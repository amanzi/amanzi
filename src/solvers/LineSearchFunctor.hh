/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//
// A functor used in line search algorithms.
//
#pragma once

namespace Amanzi {
namespace AmanziSolvers {

// functor for minimization
template<class Vector>
struct LineSearchFunctor {

  LineSearchFunctor(const Teuchos::RCP<SolverFnBase<Vector>>& my_fn,
                    const Teuchos::RCP<const Vector>& u0_,
                    const Teuchos::RCP<const Vector>& du_,
                    const Teuchos::RCP<Vector>& u_,
                    const Teuchos::RCP<Vector>& r_)
    : fn(my_fn),
      u(u_),
      u0(u0_),
      du(du_),
      r(r_),
      error(-1),
      fun_calls(0)
    {}

  double operator()(double x) const {
    u->Update(-x, *du, 1., *u0, 0);
    fn->ChangedSolution();
    fn->Residual(u, r);
    error = fn->ErrorNorm(u, r);
    fun_calls++;
    return error;
  }

  mutable double error;
  mutable int fun_calls;
  Teuchos::RCP<const Vector> u0, du;
  Teuchos::RCP<Vector> u, r;
  Teuchos::RCP<SolverFnBase<Vector>> fn;
};


} // namespace AmanziSolvers
} // namespace Amanzi
