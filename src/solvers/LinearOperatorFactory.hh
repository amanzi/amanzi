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

#include "errors.hh"
#include "exceptions.hh"

#include "LinearOperator.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "LinearOperatorNKA.hh"

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
      Teuchos::RCP<const Matrix> m,  // Apply() is required
      Teuchos::RCP<const Matrix> h);  // ApplyInverse() is required

  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> > Create(
      const string& name,
      const Teuchos::ParameterList& solvers_list,
      Teuchos::RCP<const Matrix> m) {
    return Create(name, solvers_list, m, m);
  }

  // pull name from single solver list
  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> > Create(
      Teuchos::ParameterList& slist,
      Teuchos::RCP<const Matrix> m,
      Teuchos::RCP<const Matrix> h);

  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> > Create(
      Teuchos::ParameterList& slist,
      Teuchos::RCP<const Matrix> m) {
    return Create(slist, m, m);
  }

};


/* ******************************************************************
 * The following calls have to be supported: m->Apply(...) and
 * h->ApplyInverse(...).
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> >
LinearOperatorFactory<Matrix, Vector, VectorSpace>::Create(
    const string& name, const Teuchos::ParameterList& solvers_list,
    Teuchos::RCP<const Matrix> m,
    Teuchos::RCP<const Matrix> h)
{
  if (solvers_list.isSublist(name)) {
    Teuchos::ParameterList slist = solvers_list.sublist(name);
    return Create(slist, m, h);
  } else {
    std::stringstream msgstream;
    msgstream << "LinearOperatorFactory: solver \"" << name << "\" is not on the list of solvers.";
    Errors::Message msg(msgstream.str());
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * The following calls have to be supported: m->Apply(...) and
 * h->ApplyInverse(...).
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace> >
LinearOperatorFactory<Matrix, Vector, VectorSpace>::Create(
    Teuchos::ParameterList& slist,
    Teuchos::RCP<const Matrix> m,
    Teuchos::RCP<const Matrix> h)
{
  if (slist.isParameter("iterative method")) {
    std::string method_name = slist.get<string>("iterative method");

    if (method_name == "pcg") {
      Teuchos::RCP<LinearOperatorPCG<Matrix, Vector, VectorSpace> >
          lin_op = Teuchos::rcp(new LinearOperatorPCG<Matrix, Vector, VectorSpace>(m, h));
      lin_op->Init(slist);
      lin_op->set_name(method_name);
      return lin_op;
    } else if (method_name == "gmres") {
      Teuchos::RCP<LinearOperatorGMRES<Matrix, Vector, VectorSpace> >
          lin_op = Teuchos::rcp(new LinearOperatorGMRES<Matrix, Vector, VectorSpace>(m, h));
      lin_op->Init(slist);
      lin_op->set_name(method_name);
      return lin_op;
    } else if (method_name == "nka") {
      Teuchos::RCP<LinearOperatorNKA<Matrix, Vector, VectorSpace> >
          lin_op = Teuchos::rcp(new LinearOperatorNKA<Matrix, Vector, VectorSpace>(m, h));
      lin_op->Init(slist);
      lin_op->set_name(method_name);
      return lin_op;
    } else {
      Errors::Message msg("LinearOperatorFactory: wrong value of parameter `\"iterative method`\"");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("LinearOperatorFactory: parameter `\"iterative method`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
