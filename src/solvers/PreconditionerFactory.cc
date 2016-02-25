/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base factory for preconditioners.
*/

#include "Teuchos_RCP.hpp"

#include "errors.hh"

#include "Preconditioner.hh"
#include "PreconditionerBlockILU.hh"
#include "PreconditionerBoomerAMG.hh"
#include "PreconditionerEuclid.hh"
#include "PreconditionerDiagonal.hh"
#include "PreconditionerFactory.hh"
#include "PreconditionerIdentity.hh"
#include "PreconditionerML.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
Teuchos::RCP<Preconditioner> PreconditionerFactory::Create(
    const std::string& name, const Teuchos::ParameterList& prec_list)
{
  if (prec_list.isSublist(name)) {
    Teuchos::ParameterList slist = prec_list.sublist(name);
    return Create(slist);
  } else {
    Teuchos::RCP<PreconditionerIdentity> prec = Teuchos::rcp(new PreconditionerIdentity());
    prec->Init(name, prec_list);
    return prec;
  }
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
Teuchos::RCP<Preconditioner> 
PreconditionerFactory::Create(Teuchos::ParameterList& slist)
{
  if (slist.isParameter("preconditioner type")) {
    std::string type = slist.get<std::string>("preconditioner type");

    if (type == "boomer amg") {
      Teuchos::ParameterList hypre_list = slist.sublist("boomer amg parameters");
      Teuchos::RCP<PreconditionerBoomerAMG> prec = Teuchos::rcp(new PreconditionerBoomerAMG());
      prec->Init(type, hypre_list);
      return prec;
    } else if (type == "euclid") {
      Teuchos::ParameterList hypre_list = slist.sublist("euclid parameters");
      Teuchos::RCP<PreconditionerEuclid> prec = Teuchos::rcp(new PreconditionerEuclid());
      prec->Init(type, hypre_list);
      return prec;
    } else if (type == "ml") {
      Teuchos::ParameterList ml_list = slist.sublist("ml parameters");
      Teuchos::RCP<PreconditionerML> prec = Teuchos::rcp(new PreconditionerML());
      prec->Init(type, ml_list);
      return prec;
    } else if (type == "block ilu") {
      Teuchos::ParameterList ilu_list = slist.sublist("block ilu parameters");
      Teuchos::RCP<PreconditionerBlockILU> prec = Teuchos::rcp(new PreconditionerBlockILU());
      prec->Init(type, ilu_list);
      return prec;
    } else if (type == "diagonal") {
      Teuchos::RCP<PreconditionerDiagonal> prec = Teuchos::rcp(new PreconditionerDiagonal());
      prec->Init(type, slist);
      return prec;
    } else if (type == "identity") {  // Identity preconditioner is default.
      Teuchos::RCP<PreconditionerIdentity> prec = Teuchos::rcp(new PreconditionerIdentity());
      prec->Init(type, slist);
      return prec;
    } else {
      Errors::Message msg("PreconditionerFactory: wrong value of parameter `\"preconditioner type`\"");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("PreconditionerFactory: parameter `\"preconditioner type`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}


}  // namespace AmanziPreconditioners
}  // namespace Amanzi
