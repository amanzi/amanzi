/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon
*/

#include "PK_Default.hh"

namespace Amanzi {

// Required constructor for use by the PK factory.
PK_Default::PK_Default(const Comm_ptr_type& comm,
                       Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                       const Teuchos::RCP<State>& S)
  : comm_(comm),
    name_(Keys::cleanPListName(pk_tree)),
    tag_current_(Tags::CURRENT),
    tag_next_(Tags::NEXT),
    S_(S)
{
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_plist, "PKs", true);
  if (pk_list->isSublist(name_)) {
    plist_ = Teuchos::sublist(pk_list, name_);
  } else {
    Errors::Message msg;
    msg << "There is no sublist for PK \"" << name_ << "\" in PKs list";
      Exceptions::amanzi_throw(msg);
  }

  // construct the VerboseObject
  // Note, this allows for overriding the vo plist for individual PKs in a
  // collection of PKs
  Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
  if (plist_->isSublist(name_ + " verbose object")) {
    vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
    vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
  }

  //  some tests provide nullptr
  vo_ = Teuchos::rcp(new VerboseObject(comm_, name_, *vo_plist));
};

} // namespace Amanzi
