/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Factory of linear operators.
Usage: Create("GMRES with Hypre AMG", solvers_list);
*/

#ifndef AMANZI_LINEAR_OPERATOR_FACTORY_HH_
#define AMANZI_LINEAR_OPERATOR_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"

#include "LinearOperator.hh"
#include "PCG_Operator.hh"
#include "GMRES_Operator.hh"
 
namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorFactory {
 public:
  LinearOperatorFactory() {};
  ~LinearOperatorFactory() {};

  // use name in the solvers list to initialize the solver
  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> > Create(
      const string& name, 
      const Teuchos::ParameterList& solvers_list,
      Teuchos::RCP<const Matrix> m);
};


template<class Matrix, class Vector, class VectorSpace>
Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> > 
    LinearOperatorFactory<Matrix, Vector, VectorSpace>::Create(
    const string& name, const Teuchos::ParameterList& solvers_list,
    Teuchos::RCP<const Matrix> m)
{
  if (solvers_list.isSublist(name)) {
    Teuchos::ParameterList slist = solvers_list.sublist(name);
    if (slist.isParameter("iterative method")) {
      std::string method_name = slist.get<string>("iterative method");

      if (method_name == "pcg") {
         Teuchos::RCP<PCG_Operator<Matrix, Vector, VectorSpace> > 
             lin_op = Teuchos::rcp(new PCG_Operator<Matrix, Vector, VectorSpace>(m));
         lin_op->Init(slist);
         lin_op->name() = method_name;
         return lin_op;
      } else if (method_name == "gmres") {
         Teuchos::RCP<GMRES_Operator<Matrix, Vector, VectorSpace> > 
             lin_op = Teuchos::rcp(new GMRES_Operator<Matrix, Vector, VectorSpace>(m));
         lin_op->Init(slist);
         lin_op->name() = method_name;
         return lin_op;
      } else {
        Errors::Message msg("LinearOperatorFactory: wrong value of parameter `\"iterative method`\"");
        Exceptions::amanzi_throw(msg);
      }
    } else {
      Errors::Message msg("LinearOperatorFactory: parameter `\"iterative method`\" is missing");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("LinearOperatorFactory: solver is not on the list of solvers.");
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi

#endif


