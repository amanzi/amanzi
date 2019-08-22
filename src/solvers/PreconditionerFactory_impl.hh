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

#include "PreconditionerIdentity.hh"
#include "PreconditionerDiagonal.hh"

#ifdef HAVE_EPETRA_PRECONDITIONERS
#ifdef HAVE_TRILINOS_PRECONDITIONERS
#include "PreconditionerBlockILU.hh"
#include "PreconditionerML.hh"
#endif

#ifdef HAVE_HYPRE_PRECONDITIONERS
#include "PreconditionerBoomerAMG.hh"
#include "PreconditionerEuclid.hh"
#endif
#endif


namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
template<class Matrix, class Vector>
Teuchos::RCP<Preconditioner<Matrix,Vector> >
PreconditionerFactory<Matrix,Vector>::Create(
    const std::string& name, const Teuchos::ParameterList& prec_list)
{
  if (prec_list.isSublist(name)) {
    Teuchos::ParameterList slist = prec_list.sublist(name);
    return Create(slist);
  } else {
    auto prec = Teuchos::rcp(new PreconditionerIdentity<Matrix,Vector>());
    prec->Init(name, prec_list);
    return prec;
  }
}


#ifdef HAVE_EPETRA_PRECONDITIONERS

/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
template<>
Teuchos::RCP<Preconditioner<Epetra_RowMatrix,Epetra_MultiVector> >
PreconditionerFactory<Epetra_RowMatrix,Epetra_MultiVector>::Create(Teuchos::ParameterList& slist)
{
  if (slist.isParameter("preconditioner type")) {
    std::string type = slist.get<std::string>("preconditioner type");
    
    if (type == "identity") {  // Identity preconditioner is default.
      auto prec = Teuchos::rcp(new PreconditionerIdentity<Epetra_RowMatrix,Epetra_MultiVector>());
      prec->Init(type, slist);
      return prec;

    } else if (type == "diagonal") {
      auto prec = Teuchos::rcp(new PreconditionerDiagonal<Epetra_RowMatrix,Epetra_MultiVector>());
      prec->Init(type, slist);
      return prec;

    } else if (type == "boomer amg") {
#ifdef HAVE_HYPRE_PRECONDITIONERS
      Teuchos::ParameterList hypre_list = slist.sublist("boomer amg parameters");
      auto prec = Teuchos::rcp(new PreconditionerBoomerAMG());
      prec->Init(type, hypre_list);
      return prec;
#else
      Errors::Message msg("Hypre (BoomerAMG) is not available in this installation of Amanzi.  To use Hypre, please reconfigure.");
      Exceptions::amanzi_throw(msg);
#endif

    } else if (type == "euclid") {
#ifdef HAVE_HYPRE_PRECONDITIONERS
      Teuchos::ParameterList hypre_list = slist.sublist("euclid parameters");
      auto prec = Teuchos::rcp(new PreconditionerEuclid());
      prec->Init(type, hypre_list);
      return prec;
#else
      Errors::Message msg("Hypre (Euclid) is not available in this installation of Amanzi.  To use Hypre, please reconfigure.");
      Exceptions::amanzi_throw(msg);
#endif

    } else if (type == "ml") {
#ifdef HAVE_TRILINOS_PRECONDITIONERS      
      Teuchos::ParameterList ml_list = slist.sublist("ml parameters");
      auto prec = Teuchos::rcp(new PreconditionerML());
      prec->Init(type, ml_list);
      return prec;
#else
      Errors::Message msg("ML is not available in this installation of Amanzi, this is developer error.");
      Exceptions::amanzi_throw(msg);
#endif

    } else if (type == "block ilu") {
#ifdef HAVE_TRILINOS_PRECONDITIONERS      
      Teuchos::ParameterList ilu_list = slist.sublist("block ilu parameters");
      auto prec = Teuchos::rcp(new PreconditionerBlockILU());
      prec->Init(type, ilu_list);
      return prec;
#else
      Errors::Message msg("ML is not available in this installation of Amanzi, this is developer error.");
      Exceptions::amanzi_throw(msg);
#endif

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


#endif // HAVE_EPETRA_PRECONDITIONERS


#ifdef HAVE_TPETRA_PRECONDITIONERS

/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
template<>
Teuchos::RCP<Preconditioner<Matrix_type, VectorWrapper<Vector_type> > >
PreconditionerFactory<Matrix_type,VectorWrapper<Vector_type> >::Create(Teuchos::ParameterList& slist)
{
  if (slist.isParameter("preconditioner type")) {
    std::string type = slist.get<std::string>("preconditioner type");

    if (type == "identity") {  // Identity preconditioner is default.
      auto prec = Teuchos::rcp(new PreconditionerIdentity<Matrix_type,VectorWrapper<Vector_type> >());
      prec->Init(type, slist);
      return prec;

    } else if (type == "diagonal") {
      auto prec = Teuchos::rcp(new PreconditionerDiagonal<Matrix_type,VectorWrapper<Vector_type> >());
      prec->Init(type, slist);
      return prec;

    } else if (type == "boomer amg" || type == "euclid" || type == "ml" || type == "block ilu") {
      Errors::Message msg;
      msg << "Preconditioner type \"" << type << "\" is not available in a Tpetra-based installation of Amanzi.";
      Exceptions::amanzi_throw(msg);
      
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

#endif // HAVE_TPETRA_PRECONDITIONERS

}  // namespace AmanziPreconditioners
}  // namespace Amanzi
