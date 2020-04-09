/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! The interface for a Process Kernel, an equation or system of equations.

#include "boost/algorithm/string.hpp"

#include "Key.hh"
#include "VerboseObject.hh"
#include "State.hh"

#include "PK_Default.hh"

namespace Amanzi {

PK_Default::PK_Default(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                       const Teuchos::RCP<State>& S)
  : S_(S), name_(Keys::cleanPListName(pk_tree->name()))
{
  // grab my sublist
  Teuchos::RCP<Teuchos::ParameterList> pks_list =
    Teuchos::sublist(global_plist, "PKs");
  if (pks_list->isSublist(name_)) {
    plist_ = Teuchos::sublist(pks_list, name_);
  } else {
    Errors::Message message;
    message << "PK_Default: There is no sublist for PK \"" << name_
            << "\" in PKs list\n";
    Exceptions::amanzi_throw(message);
  }

  // create a verbose object
  vo_ = Teuchos::rcp(new VerboseObject(name_, *plist_));
};

} // namespace Amanzi
