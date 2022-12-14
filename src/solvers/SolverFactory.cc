/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A factory for creating nonlinear solvers.
#include "SolverFactory.hh"

#include "SolverAA.hh"
#include "SolverNKA.hh"
#include "SolverNKA_LS.hh"
#include "SolverNKA_BT_ATS.hh"
#include "SolverNKA_LS_ATS.hh"
#include "SolverNewton.hh"
#include "SolverNox.hh"
#include "SolverJFNK.hh"
#include "SolverContinuation.hh"
#include "SolverBT.hh"

#include "Epetra_MultiVector.h"
#include "CompositeVector.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
* Initialization of the Solver
****************************************************************** */
template <class Vector, class VectorSpace>
Teuchos::RCP<Solver<Vector, VectorSpace>>
SolverFactory<Vector, VectorSpace>::Create(const std::string& name,
                                           const Teuchos::ParameterList& solver_list)
{
  if (solver_list.isSublist(name)) {
    Teuchos::ParameterList slist = solver_list.sublist(name);
    return Create(slist);
  } else {
    std::stringstream estream;
    estream << "SolverFactory: nonexistent solver sublist \"" << name << "\"";
    Errors::Message msg(estream.str());
    Exceptions::amanzi_throw(msg);
    return Teuchos::null;
  }
}


/* ******************************************************************
* Initialization of the solver
****************************************************************** */
template <class Vector, class VectorSpace>
Teuchos::RCP<Solver<Vector, VectorSpace>>
SolverFactory<Vector, VectorSpace>::Create(Teuchos::ParameterList& slist)
{
  if (slist.isParameter("solver type")) {
    std::string type = slist.get<std::string>("solver type");

    if (type == "nka") {
      if (!slist.isSublist("nka parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka parameters");
      if (!nka_list.isSublist("verbose object"))
        nka_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverNKA<Vector, VectorSpace>(nka_list));
      return solver;
    } else if (type == "aa") {
      if (!slist.isSublist("aa parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList aa_list = slist.sublist("aa parameters");
      if (!aa_list.isSublist("verbose object"))
        aa_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverAA<Vector, VectorSpace>(aa_list));
      return solver;
    } else if (type == "Newton") {
      if (!slist.isSublist("Newton parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"Newton parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList newton_list = slist.sublist("Newton parameters");
      if (!newton_list.isSublist("verbose object"))
        newton_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverNewton<Vector, VectorSpace>(newton_list));
      return solver;
    } else if (type == "nka line search") {
      if (!slist.isSublist("nka line search parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka line search parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka line search parameters");
      if (!nka_list.isSublist("verbose object"))
        nka_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverNKA_LS<Vector, VectorSpace>(nka_list));
      return solver;
    } else if (type == "nka_bt_ats") {
      if (!slist.isSublist("nka_bt_ats parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka_bt_ats parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka_bt_ats parameters");
      if (!nka_list.isSublist("verbose object"))
        nka_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverNKA_BT_ATS<Vector, VectorSpace>(nka_list));
      return solver;
    } else if (type == "nka_ls_ats") {
      if (!slist.isSublist("nka_ls_ats parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nka_ls_ats parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList nka_list = slist.sublist("nka_ls_ats parameters");
      if (!nka_list.isSublist("verbose object"))
        nka_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverNKA_LS_ATS<Vector, VectorSpace>(nka_list));
      return solver;
    } else if (type == "JFNK") {
      if (!slist.isSublist("JFNK parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"JFNK parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList jfnk_list = slist.sublist("JFNK parameters");
      if (!jfnk_list.isSublist("verbose object"))
        jfnk_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverJFNK<Vector, VectorSpace>(jfnk_list));
      return solver;
    } else if (type == "continuation") {
      if (!slist.isSublist("continuation parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"continuation parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList cont_list = slist.sublist("continuation parameters");
      if (!cont_list.isSublist("verbose object"))
        cont_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverContinuation<Vector, VectorSpace>(cont_list));
      return solver;
    } else if (type == "line search") {
      if (!slist.isSublist("line search parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"line search parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList ls_list = slist.sublist("line search parameters");
      if (!ls_list.isSublist("verbose object"))
        ls_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverBT<Vector, VectorSpace>(ls_list));
      return solver;
    } else if (type == "nox") {
      if (!slist.isSublist("nox parameters")) {
        Errors::Message msg("SolverFactory: missing sublist \"nox parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList ls_list = slist.sublist("nox parameters");
      if (!ls_list.isSublist("verbose object"))
        ls_list.set("verbose object", slist.sublist("verbose object"));
      Teuchos::RCP<Solver<Vector, VectorSpace>> solver =
        Teuchos::rcp(new SolverNox<Vector, VectorSpace>(ls_list));
      return solver;
    } else {
      Errors::Message msg("SolverFactory: wrong value of parameter `\"solver type`\"");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("SolverFactory: parameter `\"solver type`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}

template class SolverFactory<Epetra_Vector, Epetra_BlockMap>;
template class SolverFactory<CompositeVector, CompositeVectorSpace>;
template class SolverFactory<TreeVector, TreeVectorSpace>;

} // namespace AmanziSolvers
} // namespace Amanzi
