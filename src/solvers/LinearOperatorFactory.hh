/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
// #include "LinearOperatorBelosGMRES.hh"
// #include "LinearOperatorAmesos.hh"
#include "LinearOperatorNKA.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Matrix, class Vector, class VectorSpace>
class LinearOperatorFactory {
 public:
  LinearOperatorFactory(){};
  ~LinearOperatorFactory(){};

  // use name in the solvers list to initialize the solver
  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace>>
  Create(const std::string& name, const Teuchos::ParameterList& solvers_list,
         Teuchos::RCP<const Matrix> m,  // Apply() is required
         Teuchos::RCP<const Matrix> h); // ApplyInverse() is required

  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace>>
  Create(const std::string& name, const Teuchos::ParameterList& solvers_list,
         Teuchos::RCP<const Matrix> m)
  {
    return Create(name, solvers_list, m, m);
  }

  // pull name from single solver list
  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace>>
  Create(Teuchos::ParameterList& slist, Teuchos::RCP<const Matrix> m,
         Teuchos::RCP<const Matrix> h);

  Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace>>
  Create(Teuchos::ParameterList& slist, Teuchos::RCP<const Matrix> m)
  {
    return Create(slist, m, m);
  }
};


/* ******************************************************************
 * The following calls have to be supported: m->apply(...) and
 * h->applyInverse(...).
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace>>
LinearOperatorFactory<Matrix, Vector, VectorSpace>::Create(
  const std::string& name, const Teuchos::ParameterList& solvers_list,
  Teuchos::RCP<const Matrix> m, Teuchos::RCP<const Matrix> h)
{
  if (solvers_list.isSublist(name)) {
    Teuchos::ParameterList slist = solvers_list.sublist(name);
    return Create(slist, m, h);
  } else {
    std::stringstream msgstream;
    msgstream << "LinearOperatorFactory: solver \"" << name
              << "\" is not on the list of solvers.";
    Errors::Message msg(msgstream.str());
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}


/* ******************************************************************
 * The following calls have to be supported: m->apply(...) and
 * h->applyInverse(...).
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
Teuchos::RCP<LinearOperator<Matrix, Vector, VectorSpace>>
LinearOperatorFactory<Matrix, Vector, VectorSpace>::Create(
  Teuchos::ParameterList& slist, Teuchos::RCP<const Matrix> m,
  Teuchos::RCP<const Matrix> h)
{
  if (slist.isParameter("iterative method")) {
    std::string method_name = slist.get<std::string>("iterative method");

    // check for optional list of parameters
    std::string tmp(method_name);
    tmp.append(" parameters");
    if (!slist.isSublist(tmp)) {
      Teuchos::ParameterList vlist;
      Teuchos::RCP<VerboseObject> vo =
        Teuchos::rcp(new VerboseObject("Solvers::Factory", vlist));
      if (vo->getVerbLevel() >= Teuchos::VERB_LOW) {
        Teuchos::OSTab tab = vo->getOSTab();
        *vo->os() << vo->color("yellow") << "Parameter sublist \"" << tmp
                  << "\" is missing, use defaults." << vo->color() << std::endl;
      }
    }
    if (!slist.sublist(tmp).isSublist("verbose object"))
      slist.sublist(tmp).set("verbose object", slist.sublist("verbose object"));

    if (method_name == "pcg") {
      Teuchos::ParameterList pcg_list = slist.sublist("pcg parameters");
      Teuchos::RCP<LinearOperatorPCG<Matrix, Vector, VectorSpace>> lin_op =
        Teuchos::rcp(new LinearOperatorPCG<Matrix, Vector, VectorSpace>(m, h));
      lin_op->Init(pcg_list);
      lin_op->set_name(method_name);
      return lin_op;
    } else if (method_name == "gmres") {
      Teuchos::ParameterList gmres_list = slist.sublist("gmres parameters");
      Teuchos::RCP<LinearOperatorGMRES<Matrix, Vector, VectorSpace>> lin_op =
        Teuchos::rcp(
          new LinearOperatorGMRES<Matrix, Vector, VectorSpace>(m, h));
      lin_op->Init(gmres_list);
      lin_op->set_name(method_name);
      return lin_op;
      // } else if (method_name == "belos gmres") {
      //   Teuchos::ParameterList gmres_list = slist.sublist("belos gmres
      //   parameters"); Teuchos::RCP<LinearOperatorBelosGMRES<Matrix, Vector,
      //   VectorSpace> >
      //       lin_op = Teuchos::rcp(new LinearOperatorBelosGMRES<Matrix,
      //       Vector, VectorSpace>(m, h));
      //   lin_op->Init(gmres_list);
      //   lin_op->set_name(method_name);
      //   return lin_op;
    } else if (method_name == "nka") {
      Teuchos::ParameterList nka_list = slist.sublist("nka parameters");
      Teuchos::RCP<LinearOperatorNKA<Matrix, Vector, VectorSpace>> lin_op =
        Teuchos::rcp(new LinearOperatorNKA<Matrix, Vector, VectorSpace>(m, h));
      lin_op->Init(nka_list);
      lin_op->set_name(method_name);
      return lin_op;
    } else {
      Errors::Message msg(
        "LinearOperatorFactory: wrong value of parameter \"iterative method\"");
      Exceptions::amanzi_throw(msg);
    }
    // }
    // else if (slist.isParameter("direct method")) {
    //   std::string method_name = slist.get<std::string>("direct method");

    //   std::string tmp(method_name);
    //   tmp.append(" parameters");

    //   Teuchos::ParameterList amesos_list = slist.sublist(tmp);
    //   Teuchos::RCP<LinearOperatorAmesos<Matrix, Vector, VectorSpace> >
    //      lin_op = Teuchos::rcp(new LinearOperatorAmesos<Matrix, Vector,
    //      VectorSpace>(m, h));
    //   lin_op->Init(amesos_list);
    //   return lin_op;
  } else {
    Errors::Message msg("LinearOperatorFactory: parameter \"iterative method\" "
                        "or \"direct method\" not found");
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
