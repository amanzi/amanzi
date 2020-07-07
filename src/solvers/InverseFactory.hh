/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Factory of linear operators.
*/

#ifndef AMANZI_LINEAR_OPERATOR_FACTORY_HH_
#define AMANZI_LINEAR_OPERATOR_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "exceptions.hh"

#include "Inverse.hh"
#include "InverseAssembled.hh"
#include "IterativeMethodPCG.hh"
#include "IterativeMethodGMRES.hh"
//#include "InverseBelosGMRES.hh"
#include "DirectMethodAmesos.hh"
#include "DirectMethodAmesos2.hh"
#include "IterativeMethodNKA.hh"

#include "PreconditionerIdentity.hh"
#include "PreconditionerDiagonal.hh"
#include "PreconditionerIfpack.hh"

namespace Amanzi {
namespace AmanziSolvers {

//
// Helper function
//
inline Teuchos::ParameterList&
getMethodSublist(Teuchos::ParameterList& inv_list,
        const std::string& method_name)
{
  std::string method_name_pars = method_name + " parameters";

  // check for optional list of parameters, warn if not found
  if (!inv_list.isSublist(method_name_pars)) {
    Teuchos::ParameterList vlist;
    Teuchos::RCP<VerboseObject> vo = Teuchos::rcp(new VerboseObject("Solvers::Factory", vlist));
    if (vo->os_OK(Teuchos::VERB_LOW)) {
      Teuchos::OSTab tab = vo->getOSTab();
      *vo->os() << vo->color("yellow") << "Parameter sublist \"" << method_name_pars 
                << "\" is missing, using defaults." << vo->reset() << std::endl;
    }
  }
  auto& method_list = inv_list.sublist(method_name_pars);

  // update the list with verbose object info from the parent list
  if (!method_list.isSublist("verbose object"))
    method_list.set("verbose object", inv_list.sublist("verbose object"));

  // update the list with the method name
  if (!method_list.isParameter("method"))
    method_list.set<std::string>("method", method_name);
  return method_list;  
}


//
// Low level factory -- just make the object
// -----------------------------------------------------------------------------

//
// Low level factory for Iterative methods
//
template<class Matrix,
         class Preconditioner,
         class Vector=typename Matrix::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>         
Teuchos::RCP<Inverse<Matrix,Preconditioner,Vector,VectorSpace>>
createIterativeMethod(Teuchos::ParameterList& inv_list)
{
  if (!inv_list.isParameter("iterative method")) {
    Errors::Message msg("InverseFactory: parameter `\"iterative method`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  std::string method_name = inv_list.get<std::string>("iterative method");
  auto& method_list = getMethodSublist(inv_list, method_name);

  Teuchos::RCP<Inverse<Matrix,Preconditioner,Vector,VectorSpace>> inv;
  if (method_name == "gmres") {
    inv = Teuchos::rcp(new IterativeMethodGMRES<Matrix,Preconditioner,Vector,VectorSpace>());
  } else if (method_name == "pcg") {
    inv = Teuchos::rcp(new IterativeMethodPCG<Matrix,Preconditioner,Vector,VectorSpace>());
  } else if (method_name == "nka") {
    inv = Teuchos::rcp(new IterativeMethodNKA<Matrix,Preconditioner,Vector,VectorSpace>());
  // } else if (method_name == "belos gmres") {
  //   inv = Teuchos::rcp(new InverseBelosGMRES<Matrix,Preconditioner,Vector,VectorSpace>());
  } else {
    Errors::Message msg;
    msg << "Iterative method \"" << method_name << "\" is not a valid name.";
    Exceptions::amanzi_throw(msg);
  }

  if (inv.get()) inv->InitInverse(method_list);
  return inv;
}


//
// Low level factory for Direct methods on matrices that can self-assemble.
//
template<class Matrix,
         class Preconditioner,
         class Vector=typename Matrix::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>         
Teuchos::RCP<Inverse<Matrix,Preconditioner,Vector,VectorSpace>>
createDirectMethod(Teuchos::ParameterList& inv_list)
{
  if (!inv_list.isParameter("direct method")) {
    Errors::Message msg("InverseFactory: parameter `\"direct method`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  std::string method_name = inv_list.get<std::string>("direct method");

  // -- currently all direct methods require an assembled matrix
  // -- direct methods are parsed first, so we can just use list as-is
  auto inv = Teuchos::rcp(new InverseAssembled<Matrix,Preconditioner,Vector,VectorSpace>());
  inv->InitInverse(inv_list);
  return inv;
}


//
// Low level factory for Direct methods on Epetra_CrsMatrix
//
template<>
Teuchos::RCP<Inverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>>
inline
createDirectMethod(Teuchos::ParameterList& inv_list)
{
  if (!inv_list.isParameter("direct method")) {
    Errors::Message msg("InverseFactory: parameter `\"direct method`\" is missing");
    Exceptions::amanzi_throw(msg);
  }

  std::string method_name = inv_list.get<std::string>("direct method");
  auto& method_list = getMethodSublist(inv_list, method_name);

  Teuchos::RCP<Inverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>> inv = Teuchos::null;
  if (Keys::starts_with(method_name, "amesos")) {
    inv = Teuchos::rcp(new DirectMethodAmesos());
  } else if (Keys::starts_with(method_name, "amesos2")) {
    inv = Teuchos::rcp(new DirectMethodAmesos2());
  } else {
    Errors::Message msg;
    msg << "Direct method \"" << method_name << "\" is not a valid name. Currently only \"amesos: *\" or \"amesos2: *\" are valid options.";
    Exceptions::amanzi_throw(msg);
  }
  if (inv.get()) inv->InitInverse(method_list);
  return inv;
}


//
// Low level factory for preconditioning methods on matrices that can self-assemble
//
template<class Matrix,
         class Preconditioner,
         class Vector=typename Matrix::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>         
Teuchos::RCP<Inverse<Matrix,Preconditioner,Vector,VectorSpace>>
createPreconditioner(Teuchos::ParameterList& inv_list)
{
  if (!inv_list.isParameter("preconditioning method")) {
    Errors::Message msg("PreconditionerFactory: parameter `\"preconditioning method`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  std::string method_name = inv_list.get<std::string>("preconditioning method");

  Teuchos::RCP<Inverse<Matrix,Preconditioner,Vector,VectorSpace>> inv = Teuchos::null;
  if (method_name == "identity") {
    // NOTE: add diagonal to this once operators have been refactored to
    // provide their diagonal without assembly, see #455
    auto& method_list = getMethodSublist(inv_list, method_name);
    inv = Teuchos::rcp(new PreconditionerIdentity<Matrix,Preconditioner,Vector,VectorSpace>());
    inv->InitInverse(method_list);
  } else {
    // other preconditioners require assembly
    inv = Teuchos::rcp(new InverseAssembled<Matrix,Preconditioner,Vector,VectorSpace>());
    inv->InitInverse(inv_list);
  }
  return inv;
}


//
// Low level factory for preconditioning methods on Epetra_CrsMatrix
//
template<>
Teuchos::RCP<Inverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>>
inline
createPreconditioner(Teuchos::ParameterList& inv_list)
{
  if (!inv_list.isParameter("preconditioning method")) {
    Errors::Message msg("PreconditionerFactory: parameter `\"preconditioning method`\" is missing");
    Exceptions::amanzi_throw(msg);
  }

  std::string method_name = inv_list.get<std::string>("preconditioning method");
  auto& method_list = getMethodSublist(inv_list, method_name);

  Teuchos::RCP<Inverse<Epetra_CrsMatrix,Epetra_CrsMatrix,Epetra_Vector,Epetra_Map>> inv = Teuchos::null;
  if (method_name == "diagonal") {
    inv = Teuchos::rcp(new PreconditionerDiagonal());
  } else if (method_name == "block ilu") {
    method_list.set<std::string>("method", "ILU");
    inv = Teuchos::rcp(new PreconditionerIfpack());
  } else if (Keys::starts_with(method_name, "ifpack: ")) {
    method_list.set<std::string>("method", method_name.substr(std::string("ifpack: ").length(),
            method_name.length()));
    inv = Teuchos::rcp(new PreconditionerIfpack());
    //  } else if (method_name == "boomer amg") {
  //   inv = Teuchos::rcp(new PreconditionerBoomerAMG());
  // } else if (method_name == "euclid") {
  //   inv = Teuchos::rcp(new PreconditionerEuclid());
  // } else if (method_name == "ml") {
  //   inv = Teuchos::rcp(new PreconditionerML());
  // } else if (method_name == "block ilu") {
  //   inv = Teuchos::rcp(new PreconditionerBlockILU());
  // } else if (method_name == "diagonal") {
  //   inv = Teuchos::rcp(new PreconditionerAssembledDiagonal());
  } else if (method_name == "identity") {
    inv = Teuchos::rcp(new PreconditionerIdentity<Epetra_CrsMatrix,Epetra_CrsMatrix,
                       Epetra_Vector,Epetra_Map>());
  } else {
    Errors::Message msg;
    msg << "Preconditioning method \"" << method_name << "\" is not a valid name.";
    Exceptions::amanzi_throw(msg);
  }
  if (inv.get()) inv->InitInverse(method_list);
  return inv;
}


//
// Higher level method that parses all of iterative, direct, and precondioners,
// and generates a single object that does both forward and inverse apply.
//
// NOTE: the return type of this is the more general one!
template<class Operator,
         class Preconditioner,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(Teuchos::ParameterList& inv_list,
              const Teuchos::RCP<Operator>& m,
              const Teuchos::RCP<Preconditioner>& h)
{
  Teuchos::RCP<Matrix<Vector,VectorSpace>> inv;
  if (inv_list.isParameter("direct method")) {
    auto dir_inv = createDirectMethod<Operator,Preconditioner,Vector,VectorSpace>(inv_list);
    dir_inv->set_matrices(m,h);
    inv = dir_inv;
    return inv;
  }

  if (inv_list.isParameter("preconditioning method")) {
    auto pc_inv = createPreconditioner<Operator,Preconditioner,Vector,VectorSpace>(inv_list);
    pc_inv->set_matrices(m,h);
    inv = pc_inv;
  
    if (inv_list.isParameter("iterative method")) {
      auto inv_iter = createIterativeMethod<Operator,typeof(*inv),Vector,VectorSpace>(inv_list);
      inv_iter->set_matrices(m, inv);
      inv = inv_iter;
    }
    
  } else {
    if (inv_list.isParameter("iterative method")) {
      auto inv_iter = createIterativeMethod<Operator,Preconditioner,Vector,VectorSpace>(inv_list);
      inv_iter->set_matrices(m, h);
      inv = inv_iter;
    }
  }
  return inv;
}


//
// Factory style preferred in Amanzi -- a global, const plist of options,
// indexed by a local name.
//
template<class Operator,
         class Preconditioner,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(const std::string& name,
              const Teuchos::ParameterList& solvers_list,
              const Teuchos::RCP<Operator>& m,
              const Teuchos::RCP<Preconditioner>& h)
{
  if (solvers_list.isSublist(name)) {
    Teuchos::ParameterList inv_list = solvers_list.sublist(name);
    return createInverse<Operator,Preconditioner,Vector,VectorSpace>(
        inv_list, m, h);
  } else {
    Errors::Message msg;
    msg << "InverseFactory: method \"" << name << "\" is not available.";
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}

template<class Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>         
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(const std::string& name,
              const Teuchos::ParameterList& solvers_list,
              const Teuchos::RCP<Operator>& m)
{
  return createInverse<Operator,Operator,Vector,VectorSpace>(
      name, solvers_list, m, m);
}


//
// Factory style preferred in ATS -- a local list modified in place to preserve
// changes.
//
template<class Operator,
         class Preconditioner,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>         
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(Teuchos::ParameterList& inv_list,
              const Teuchos::RCP<const Operator>& m,
              const Teuchos::RCP<const Preconditioner>& h)
{
  auto inv = createInverse<Operator,Preconditioner,Vector,VectorSpace>(inv_list);
  inv->set_matrices(m,h);
  return inv;
}  
    
template<class Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>         
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(Teuchos::ParameterList& inv_list,
              const Teuchos::RCP<const Operator>& m)
{
  return createInverse<Operator,Operator,Vector,VectorSpace>(inv_list, m, m);
}


}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
