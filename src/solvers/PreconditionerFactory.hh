/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_PRECONDITIONER_FACTORY_HH_
#define AMANZI_PRECONDITIONER_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Preconditioner.hh"
#include "PreconditionerIdentity.hh"
#include "PreconditionerDiagonal.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template <class Matrix, class Vector>
class PreconditionerFactory {
 public:
  PreconditionerFactory(){};
  ~PreconditionerFactory(){};

  Teuchos::RCP<Preconditioner<Matrix, Vector>>
  Create(const std::string& name, const Teuchos::ParameterList& prec_list);
  Teuchos::RCP<Preconditioner<Matrix, Vector>>
  Create(Teuchos::ParameterList& prec_list);
};


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
template <class Matrix, class Vector>
Teuchos::RCP<Preconditioner<Matrix, Vector>>
PreconditionerFactory<Matrix, Vector>::Create(
  const std::string& name, const Teuchos::ParameterList& prec_list)
{
  if (prec_list.isSublist(name)) {
    Teuchos::ParameterList slist = prec_list.sublist(name);
    return Create(slist);
  } else {
    auto prec = Teuchos::rcp(new PreconditionerIdentity<Matrix, Vector>());
    prec->Init(name, prec_list);
    return prec;
  }
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
template <class Matrix, class Vector>
Teuchos::RCP<Preconditioner<Matrix, Vector>>
PreconditionerFactory<Matrix, Vector>::Create(Teuchos::ParameterList& slist)
{
  if (slist.isParameter("preconditioner type")) {
    std::string type = slist.get<std::string>("preconditioner type");

    if (type == "identity") { // Identity preconditioner is default.
      auto prec = Teuchos::rcp(new PreconditionerIdentity<Matrix, Vector>());
      prec->Init(type, slist);
      return prec;

    } else if (type == "diagonal") {
      auto prec = Teuchos::rcp(new PreconditionerDiagonal<Matrix, Vector>());
      prec->Init(type, slist);
      return prec;

    } else if (type == "boomer amg" || type == "euclid" || type == "ml" ||
               type == "block ilu") {
      Errors::Message msg;
      msg << "Preconditioner type \"" << type
          << "\" is not available in a Tpetra-based installation of Amanzi.";
      Exceptions::amanzi_throw(msg);

    } else {
      Errors::Message msg("PreconditionerFactory: wrong value of parameter "
                          "`\"preconditioner type`\"");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg(
      "PreconditionerFactory: parameter `\"preconditioner type`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}


} // namespace AmanziPreconditioners
} // namespace Amanzi


#endif
