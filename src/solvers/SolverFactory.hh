/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A factory for creating nonlinear solvers.
/*!

Nonlinear solvers are used within implicit time integration schemes to drive
the residual to zero and thereby solve for the primary variable at the new
time.

.. _solver-typed-spec:
.. admonition:: solver-typed-spec

    * `"solver type`" ``[string]`` Type of the solver.  One of:

      - `"Newton`" See `Solver: Newton and Inexact Newton`_
      - `"JFNK`" See `Solver: Jacobian-Free Newton Krylov`_
      - `"line search`" See `Solver: Newton with Line Search`_
      - `"continuation`" See `Solver: Nonlinear Continuation`_
      - `"nka`" See `Solver: Nonlinear Krylov Acceleration`_
      - `"aa`" See `Solver: Anderson Acceleration`_
      - `"nka line search`" See `Solver: NKA with Line Search`"
      - `"nka_ls_ats`" See `Solver: NKA with Line Search, ATS`_
      - `"nka_bt_ats`" See `Solver: NKA with backtracking, ATS`_
      - `"nox`" See `Solver: NOX`_

    * `"_solver_type_ parameters`" ``[_solver_type_-spec]`` A sublist containing
      parameters specific to the type.

.. warning::

    `"JFNK`", `"line search`", and `"continuation`" methods have not been
    beaten on as much as other methods.  `"nka_ls_ats`" is somewhat deprecated
    and probably shouldn't be used.  Prefer `"nka`" for simple problems,
    `"nka_bt_ats`" for freeze-thaw problems or other problems with strong
    nonlinearities, and `"Newton`" when you have a good Jacobian.  While
    `"nox`" hasn't been used extensively, it may be quite useful.



*/

#ifndef AMANZI_SOLVER_FACTORY_HH_
#define AMANZI_SOLVER_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "Solver.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverFactory {
 public:
  Teuchos::RCP<Solver<Vector, VectorSpace>>
  Create(const std::string& name, const Teuchos::ParameterList& solver_list);

  Teuchos::RCP<Solver<Vector, VectorSpace>> Create(Teuchos::ParameterList& solver_list);
};

} // namespace AmanziSolvers
} // namespace Amanzi


#endif
