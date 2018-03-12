/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Factory for nonlinear solvers.
*/

#ifndef AMANZI_SOLVER_FACTORY_HH_
#define AMANZI_SOLVER_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "Solver.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector,class VectorSpace>
struct SolverFactory {
 public:
  Teuchos::RCP<Solver<Vector,VectorSpace> >
  Create(const std::string& name,
         const Teuchos::ParameterList& solver_list);

  Teuchos::RCP<Solver<Vector,VectorSpace> >
  Create(Teuchos::ParameterList& solver_list);
};

}  // namespace Amanzisolvers
}  // namespace Amanzi


#endif
