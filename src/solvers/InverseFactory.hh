/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

/*

Factory functions for creating Inverse objects which implement the
ApplyInverse() method.


Note that we split these into three approaches --

- direct methods, e.g. LU, which expect an assembled matrix, provide a true
  inverse, and do not use a preconditioner

- iterative methods, e.g. GMRES, which do not need an assembled matrix, provide
  a true inverse, and require a preconditioner (even if just the identity)

- preconditioning methods, e.g. ILU, which sometimes but not always use an
  assembled matrix, provide an approximate inverse, and do not use a
  preconditioner themselves

From an implementation perspective, we split this factory into three cases,
depending upon the supplied Preconditioner.  Traits are provided for both
is_assembled (e.g. Epetra_CrsMatrix is itself an assembled matrix, and
therefore provides access methods) and is_assembling (e.g. Operator is
self-assembling, and therefore provides AssembleMatrix() and A() access
methods).  These traits are used to write factories that supply the correct set
of options depending upon whether the template argument is_assembled or
is_assembling.

NOTE: we are currently punting on the case that both is_assembled and
is_assembling, but this should never be true?

Finally, there are two approaches to supplying ParameterLists to this factory.
Amanzi prefers a global list of available methods, which can be indexed with a
name.  As these may be used multiple times, the global list is copied, then
passed in to allow multiple uses each getting their own copy of the list.  ATS
prefers that this copy is already in the input spec -- there is no global list
of methods.
  
*/

#ifndef AMANZI_LINEAR_OPERATOR_FACTORY_HH_
#define AMANZI_LINEAR_OPERATOR_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "exceptions.hh"

#include "InverseHelpers.hh"
#include "Inverse.hh"
#include "InverseAssembled_decl.hh"
#include "IterativeMethodPCG.hh"
#include "IterativeMethodGMRES.hh"
#include "IterativeMethodBelos.hh"
#include "DirectMethodAmesos.hh"
#include "DirectMethodAmesos2.hh"
#include "IterativeMethodNKA.hh"

#include "PreconditionerIdentity.hh"
#include "PreconditionerDiagonal.hh"
#include "PreconditionerIfpack.hh"
#include "PreconditionerHypre.hh"
#include "PreconditionerML.hh"

namespace Amanzi {
namespace AmanziSolvers {

namespace Impl {

//
// Helper functions
//

inline void
warn(const std::string& warning)
{
  Teuchos::ParameterList vlist;
  auto vo = Teuchos::rcp(new VerboseObject("Solvers::Factory", vlist));
  Teuchos::OSTab tab = vo->getOSTab();
  if (vo->os_OK(Teuchos::VERB_LOW)) {
      *vo->os() << vo->color("yellow") << warning << vo->reset() << std::endl;
  }  
}  


inline Teuchos::ParameterList&
getMethodSublist(Teuchos::ParameterList& inv_list,
        const std::string& method_name)
{
  std::string method_name_pars = method_name + " parameters";

  // check for optional list of parameters, warn if not found
  if (!inv_list.isSublist(method_name_pars)) {
    std::string warning = std::string("Parameter sublist \"")
                          + method_name_pars + "\" is missing, using defaults.";
    warn(warning);
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

} // namespace Impl


//
// Helper function that ensures that "make one iteration" is in the parameter
// list.
//
inline void
setMakeOneIterationCriteria(Teuchos::ParameterList& plist) {
  if (plist.isParameter("iterative method")) {
    auto& method_list = Impl::getMethodSublist(plist, plist.get<std::string>("iterative method"));
    Teuchos::Array<std::string> criteria;
    criteria = method_list.get<Teuchos::Array<std::string>>("convergence criteria", criteria);
    if (std::find(criteria.begin(), criteria.end(), "make one iteration") == criteria.end()) {
      criteria.push_back("make one iteration");
    }
    if (criteria.size() == 1) {
      // has only make one iteration
      criteria.push_back("relative rhs");
    }
    method_list.set("convergence criteria", criteria);
  }
}

//
// Helper function to merge Amanzi-style "solvers" and "preconditioners" lists
// options.
//
inline Teuchos::ParameterList
mergePreconditionerSolverLists(
    const std::string& pc_name, const Teuchos::ParameterList& pc_list,
    const std::string& ls_name, const Teuchos::ParameterList& ls_list,
    bool make_one_iteration=false)
{
  Teuchos::ParameterList inv_list;
  if (!ls_name.empty() && ls_name != "none") {
    if (!ls_list.isSublist(ls_name)) {
      Errors::Message msg;
      msg << "Requested solver method: \"" << ls_name
          << "\" is not a valid name provided in the \"solvers\" list.";
      Exceptions::amanzi_throw(msg);
    } else {
      inv_list.setParameters(ls_list.sublist(ls_name));
    }
  }
  if (!pc_name.empty() && pc_name != "none") {
    if (!pc_list.isSublist(pc_name)) {
      Errors::Message msg;
      msg << "Requested preconditioner method: \"" << pc_name
          << "\" is not a valid name provided in the \"preconditioners\" list.";
      Exceptions::amanzi_throw(msg);
    } else {
      inv_list.setParameters(pc_list.sublist(pc_name));
    }
  }
  if (make_one_iteration) setMakeOneIterationCriteria(inv_list);
  return inv_list;
}



//
// Low level factory -- just make the object
// -----------------------------------------------------------------------------

//
// Iterative methods work with preconditioners
//
template<class Operator,
         class Preconditioner=Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Operator::VectorSpace_t>         
Teuchos::RCP<Inverse<Operator,Preconditioner,Vector,VectorSpace>>
createIterativeMethod(const std::string& method_name,
                           Teuchos::ParameterList& inv_list)
{
  auto& method_list = Impl::getMethodSublist(inv_list, method_name);

  Teuchos::RCP<Inverse<Operator,Preconditioner,Vector,VectorSpace>> inv;
  if (method_name == "gmres") {
    inv = Teuchos::rcp(new IterativeMethodGMRES<Operator,Preconditioner,Vector,VectorSpace>());
  } else if (method_name == "pcg") {
    inv = Teuchos::rcp(new IterativeMethodPCG<Operator,Preconditioner,Vector,VectorSpace>());
  } else if (method_name == "nka") {
    inv = Teuchos::rcp(new IterativeMethodNKA<Operator,Preconditioner,Vector,VectorSpace>());
  } else if (Keys::starts_with(method_name, "belos")) {
    inv = Teuchos::rcp(new IterativeMethodBelos<Operator,Preconditioner,Vector,VectorSpace>());
  } else {
    Errors::Message msg;
    msg << "Preconditioned method \"" << method_name << "\" is not a valid name.";
    Exceptions::amanzi_throw(msg);
  }

  if (inv.get()) inv->set_inverse_parameters(method_list);
  return inv;
}

//
// This also potentially gets used by client code...
//
template<class Operator,
         class Preconditioner,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Operator::VectorSpace_t>
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createIterativeMethod(Teuchos::ParameterList& inv_list,
              const Teuchos::RCP<Operator>& m,
              const Teuchos::RCP<Preconditioner>& h)
{
  auto method_name = inv_list.get<std::string>("iterative method");
  auto inv = createIterativeMethod<Operator,Preconditioner,Vector,VectorSpace>(method_name, inv_list);
  inv->set_matrices(m, h);
  return inv;
}

template<class Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Operator::VectorSpace_t>
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createIterativeMethod(Teuchos::ParameterList& inv_list,
                      const Teuchos::RCP<Operator>& m)
{
  auto method_name = inv_list.get<std::string>("iterative method");
  auto inv = createIterativeMethod<Operator,Operator,Vector,VectorSpace>(method_name, inv_list);
  inv->set_matrices(m, m);
  return inv;
}



//
// Assembled methods work on matrices.
//
template<class Matrix=Epetra_CrsMatrix,
         class Vector=Epetra_Vector,
         class VectorSpace=Epetra_Map>
Teuchos::RCP<Inverse<Matrix,Matrix,Vector,VectorSpace>>
createAssembledMethod(const std::string& method_name, Teuchos::ParameterList& inv_list)
{
  auto& method_list = Impl::getMethodSublist(inv_list, method_name);

  Teuchos::RCP<Inverse<Matrix,Matrix,Vector,VectorSpace>> inv = Teuchos::null;
  if (Keys::starts_with(method_name, "amesos")) {
    int amesos_version = 0;
    // figure out the version
    // -- old style -- from a parameter
    if (method_list.isParameter("amesos version")) {
      amesos_version = method_list.get<int>("amesos version");
    } else if (Keys::starts_with(method_name, "amesos2")) {
      amesos_version = 2;
    } else {
      amesos_version = 1;
    }
    if (amesos_version == 1) {
      inv = Teuchos::rcp(new DirectMethodAmesos());
    } else {
      inv = Teuchos::rcp(new DirectMethodAmesos2());
    }
  } else if (method_name == "diagonal") {
    inv = Teuchos::rcp(new PreconditionerDiagonal());
  } else if (method_name == "block ilu") {
    method_list.set<std::string>("method", "ILU");
    inv = Teuchos::rcp(new PreconditionerIfpack());
  } else if (Keys::starts_with(method_name, "ifpack: ")) {
    method_list.set<std::string>("method", method_name.substr(std::string("ifpack: ").length(),
            method_name.length()));
    inv = Teuchos::rcp(new PreconditionerIfpack());
  } else if (method_name == "boomer amg" || method_name == "euclid") {
    method_list.set<std::string>("method", method_name);
    inv = Teuchos::rcp(new PreconditionerHypre());
  } else if (Keys::starts_with(method_name, "hypre: ")) {
    method_list.set<std::string>("method", method_name.substr(std::string("hypre: ").length(),
            method_name.length()));
    inv = Teuchos::rcp(new PreconditionerHypre());
  } else if (method_name == "ml") {
    inv = Teuchos::rcp(new PreconditionerML());
  } else if (method_name == "identity") {
    inv = Teuchos::rcp(new PreconditionerIdentity<Matrix,Matrix,Vector,VectorSpace>());

  } else {
    Errors::Message msg;
    msg << "Direct method \"" << method_name << "\" is not a valid name. Currently only \"amesos: *\" or \"amesos2: *\" are valid options.";
    Exceptions::amanzi_throw(msg);
  }
  if (inv.get()) inv->set_inverse_parameters(method_list);
  return inv;
}


//
// Higher level method that parses all of iterative, direct, and precondioners,
// and generates a single object that does both forward and inverse apply.
//
//
// Factory style preferred in ATS -- a local list modified in place to preserve
// changes.
//
// NOTE: the return type of this is the more general one!
template<class Operator,
         class Assembler=Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Operator::VectorSpace_t>
typename std::enable_if<Impl::is_assembling<Assembler>::value,
                        Teuchos::RCP<Matrix<Vector,VectorSpace>>>::type
createInverse(Teuchos::ParameterList& inv_list,
              const Teuchos::RCP<Operator>& m,
              const Teuchos::RCP<Assembler>& h)
{
  // deal with deprecated option
  if (inv_list.isParameter("preconditioner type")) {
    Impl::warn("InverseFactory: DEPRECATION -- please rename \"preconditioner type\" to \"preconditioning method\".  \"preconditioner type\" may silently be ignored in the future.");
    if (!inv_list.isParameter("preconditioning method"))
      inv_list.set<std::string>("preconditioning method", inv_list.get<std::string>("preconditioner type"));
  }
  
  using Matrix_t = Matrix<Vector,VectorSpace>;
  using Inverse_t = Inverse<Operator,Assembler,Vector,VectorSpace>;

  Teuchos::RCP<Matrix_t> inv = Teuchos::null;

  std::string method_name;
  if (inv_list.isParameter("direct method")) {
    method_name = inv_list.get<std::string>("direct method");
  } else if (inv_list.isParameter("preconditioning method")) {
    method_name = inv_list.get<std::string>("preconditioning method");
  } else if (inv_list.isParameter("iterative method")) {
    method_name = "identity";
    Impl::warn(std::string("InverseFactory: WARNING -- iterative method \"")+
         inv_list.get<std::string>("iterative method")+
         "\" requires a preconditioner, but no parameter \"preconditioning method\" was supplied: using \"identity\"");
  } else {
    Errors::Message msg("InverseFactory (Assembler): at least one of \"direct method\", \"iterative method\", or \"preconditioning method\" must be supplied.");
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::RCP<Inverse_t> dir_inv;
  // only "identity" is a valid "direct method" that is not assembled
  if (method_name == "identity") {
    dir_inv = Teuchos::rcp(new PreconditionerIdentity<Operator,Assembler,Vector,VectorSpace>());
  } else {
    dir_inv = Teuchos::rcp(new InverseAssembled<Operator,Assembler,Vector,VectorSpace>(method_name));
  }
  dir_inv->set_inverse_parameters(inv_list);
  dir_inv->set_matrices(m,h);
  inv = dir_inv;

  if (inv_list.isParameter("iterative method")) {
    if (inv_list.isParameter("direct method")) {
      Impl::warn("InverseFactory: WARNING -- both \"direct method\" and \"iterative method\" were supplied -- using the direct method.");
    } else {
      auto iter_inv = createIterativeMethod<Operator,Inverse_t,Vector,VectorSpace>(inv_list, m, dir_inv);
      inv = iter_inv;
    }
  }
  return inv;
}


template<class Operator,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Operator::VectorSpace_t>         
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(Teuchos::ParameterList& inv_list,
              const Teuchos::RCP<Operator>& m)
{
  return createInverse<Operator,Operator,Vector,VectorSpace>(inv_list, m, m);
}


//
// Factory style preferred in Amanzi -- a global, const plist of options,
// indexed by a local name.
//
template<class Operator,
         class Preconditioner,
         class Vector=typename Operator::Vector_t,
         class VectorSpace=typename Operator::VectorSpace_t>
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
         class VectorSpace=typename Operator::VectorSpace_t>         
Teuchos::RCP<Matrix<Vector,VectorSpace>>
createInverse(const std::string& name,
              const Teuchos::ParameterList& solvers_list,
              const Teuchos::RCP<Operator>& m)
{
  return createInverse<Operator,Operator,Vector,VectorSpace>(
      name, solvers_list, m, m);
}


}  // namespace AmanziSolvers
}  // namespace Amanzi

#include "InverseAssembled_impl.hh"

#endif
