/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Process Kernels

*/

#include "EvaluatorIndependentFunction.hh"
#include "State.hh"
#include "TreeVector.hh"
#include "PK_Helpers.hh"

#include "PK_Physical_Default.hh"

namespace Amanzi {


void
PK_Physical_Default::parseParameterList()
{
  key_ = Keys::readKey(*plist_, domain_, "primary variable");
  passwd_ = plist_->get<std::string>("primary variable password", name());

  // require primary variable evaluators
  requireEvaluatorAtNext(key_, tag_next_, *S_, passwd_);
  requireEvaluatorAtCurrent(key_, tag_current_, *S_, passwd_);
}

void
PK_Physical_Default::Setup()
{
  // set up the debugger
  Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
  if (plist_->isSublist(name_ + " verbose object")) {
    vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
    vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
  }

  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *vo_plist));
};


void
PK_Physical_Default::CommitStep(double t_old, double t_new, const Tag& tag_next)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state @ " << tag_next << std::endl;

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  assign(key_, tag_current, tag_next, *S_);
}


void
PK_Physical_Default::FailStep(double t_old, double t_new, const Tag& tag_next)
{
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  assign(key_, tag_next, tag_current, *S_);
}


// -----------------------------------------------------------------------------
//  Marks as changed
// -----------------------------------------------------------------------------
void
PK_Physical_Default::ChangedSolutionPK(const Tag& tag)
{
  changedEvaluatorPrimary(key_, tag, *S_);
}


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void
PK_Physical_Default::Initialize()
{
  // Get the record
  Record& record = S_->GetRecordW(key_, tag_next_, passwd_);

  // Initialize the data
  if (!record.initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial conditions")) {
      Errors::Message message;
      message << name() << " has no \"initial conditions\" sublist.";
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial conditions");
    record.Initialize(ic_plist);

    // communicate just to make sure values are initialized for valgrind's sake
    if (record.Get<CompositeVector>().Ghosted())
      record.Get<CompositeVector>().ScatterMasterToGhosted();
    ChangedSolutionPK(tag_next_);
  }

  if (solution_ != Teuchos::null) solution_->SetData(record.GetPtrW<CompositeVector>(passwd_));
};


} // namespace Amanzi
