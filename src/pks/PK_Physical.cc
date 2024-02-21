/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky
*/

/*
  Process Kernels

*/

#include "PK_Physical.hh"

#include "State.hh"
#include "TreeVector.hh"
#include "Debugger.hh"
#include "PK_Helpers.hh"


namespace Amanzi {

// -----------------------------------------------------------------------------
// Transfer operators -- copies ONLY pointers
// -----------------------------------------------------------------------------
Teuchos::RCP<TreeVectorSpace>
PK_Physical::getSolutionSpace() const
{
  CompositeVectorSpace& cvs = S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_);
  auto soln_space = Teuchos::rcp(new TreeVectorSpace(cvs.getComm()));
  soln_space->setData(cvs.CreateSpace());
  return soln_space;
}

void
PK_Physical::moveStateToSolution(const Tag& tag, TreeVector& solution)
{
  solution.setData(S_->GetPtrW<CompositeVector>(key_, tag, getName()));
}


void
PK_Physical::moveSolutionToState(const TreeVector& solution, const Tag& tag)
{
  AMANZI_ASSERT(solution.getData() == S_->GetPtr<CompositeVector>(key_, tag));
}


// -----------------------------------------------------------------------------
// Construction -- parse the ParameterList
// -----------------------------------------------------------------------------
void
PK_Physical::parseParameterList()
{
  domain_ = plist_->get<std::string>("domain name", "domain");
  mesh_ = S_->GetMesh(domain_);
  key_ = Keys::readKey(*plist_, domain_, "primary variable");

  PK::parseParameterList();

  // primary variable max change
  max_valid_change_ = plist_->get<double>("max valid change", -1.0);

  // set up the debugger
  Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
  if (plist_->isSublist(name_ + " verbose object")) {
    vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
    vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
  }
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *vo_plist));
}


// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------
void
PK_Physical::setup()
{
  // require primary variable evaluators
  PKHelpers::requireAtNext(key_, tag_next_, *S_, name_);
  PKHelpers::requireAtCurrent(key_, tag_current_, *S_, name_);
};


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void
PK_Physical::initialize()
{
  // Get the record
  Record& record = S_->GetRecordW(key_, tag_next_, getName());

  // Initialize the data
  if (!record.initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      Errors::Message message;
      message << getName() << " has no initial condition parameter list.";
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    ic_plist.setName(key_);
    record.Initialize(ic_plist);

    // debug
    db_->WriteVector(key_ + " IC", record.GetPtr<CompositeVector>().ptr(), true);

    // communicate just to make sure values are initialized for valgrind's sake
    record.Get<CompositeVector>().scatterMasterToGhosted();
    markChangedSolution(tag_next_);
  }
};


// -----------------------------------------------------------------------------
// Called after all PKs have finished advancing a step successfully
// -----------------------------------------------------------------------------
void
PK_Physical::commitStep(double t_old, double t_new, const Tag& tag_next)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state @ " << tag_next << std::endl;

  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(key_, tag_current, tag_next, *S_);
}


// -----------------------------------------------------------------------------
// Called if any PKs have failed to advance a step successfully
// -----------------------------------------------------------------------------
void
PK_Physical::failStep(double t_old, double t_new, const Tag& tag_next)
{
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(key_, tag_next, tag_current, *S_);
}


// -----------------------------------------------------------------------------
// Ensures the step size is smaller than max_valid_change
// -----------------------------------------------------------------------------
bool
PK_Physical::isValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Validating time step." << std::endl;

  if (max_valid_change_ > 0.0) {
    const CompositeVector& var_new = S_->Get<CompositeVector>(key_, tag_next_);
    const CompositeVector& var_old = S_->Get<CompositeVector>(key_, tag_current_);
    CompositeVector dvar(var_new);
    dvar.update(-1., var_old, 1.);
    double change = dvar.normInf();
    if (change > max_valid_change_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max primary variable change=" << change
                   << " > limit=" << max_valid_change_ << std::endl;
      return false;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
//  Marks as changed
// -----------------------------------------------------------------------------
void
PK_Physical::markChangedSolution(const Tag& tag)
{
  PKHelpers::changedEvaluatorPrimary(key_, tag, *S_);
}

} // namespace Amanzi
