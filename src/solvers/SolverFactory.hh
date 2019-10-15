/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A factory for creating nonlinear solvers.

/*!

The solver factory must be provided a list containing at least two things --
the \"solver type\" parameter and a sublist containing parameters for that
solver.

* `"solver type`" ``[string]`` Type of the solver, one of the below.
* `"_solver_type_ parameters`" ``[_solver_type_-spec]`` A sublist containing
  parameters specific to the type.

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
struct SolverFactory {
 public:
  Teuchos::RCP<Solver<Vector, VectorSpace>>
  Create(const std::string& name, const Teuchos::ParameterList& solver_list);

  Teuchos::RCP<Solver<Vector, VectorSpace>>
  Create(Teuchos::ParameterList& solver_list);
};

} // namespace AmanziSolvers
} // namespace Amanzi


#endif
