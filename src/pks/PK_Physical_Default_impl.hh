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

#pragma once

#include "PK_Physical_Default_decl.hh"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "Debugger.hh"
#include "State.hh"
#include "PK_Helpers.hh"


namespace Amanzi {

// -----------------------------------------------------------------------------
// Transfer operators -- copies ONLY pointers
// -----------------------------------------------------------------------------
template <class PK_type>
Teuchos::RCP<TreeVectorSpace>
PK_Physical_Default<PK_type>::getSolutionSpace() const
{
  CompositeVectorSpace& cvs =
    S_->template Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_);
  auto soln_space = Teuchos::rcp(new TreeVectorSpace(cvs.getComm()));
  soln_space->setData(cvs.CreateSpace());
  return soln_space;
}


template <class PK_type>
void
PK_Physical_Default<PK_type>::moveStateToSolution(const Tag& tag, TreeVector& solution)
{
  solution.setData(S_->template GetPtrW<CompositeVector>(key_, tag, PK_type::getName()));
}


template <class PK_type>
void
PK_Physical_Default<PK_type>::moveSolutionToState(const TreeVector& solution, const Tag& tag)
{
  AMANZI_ASSERT(solution.getData() == S_->template GetPtr<CompositeVector>(key_, tag));
}


// -----------------------------------------------------------------------------
// Construction -- parse the ParameterList
// -----------------------------------------------------------------------------
template <class PK_type>
void
PK_Physical_Default<PK_type>::parseParameterList()
{
  PK_type::parseParameterList();

  domain_ = plist_->template get<std::string>("domain name", "domain");
  mesh_ = S_->GetMesh(domain_);
  key_ = Keys::readKey(*plist_, domain_, "primary variable");

  // primary variable max change
  max_valid_change_ = plist_->template get<double>("max valid change", -1.0);

  // set up the debugger
  Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
  if (plist_->isSublist(name_ + " verbose object")) {
    vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
    vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
  }
  db_ = Teuchos::rcp(new Debugger(AmanziMesh::onMemHost(mesh_), name_, *vo_plist));
}


// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------
template <class PK_type>
void
PK_Physical_Default<PK_type>::setup()
{
  // require primary variable evaluators
  PKHelpers::requireAtNext(key_, tag_next_, *S_, name_);
  PKHelpers::requireAtCurrent(key_, tag_current_, *S_, name_);

  PK_type::setup();
};


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
template <class PK_type>
void
PK_Physical_Default<PK_type>::initialize()
{
  // Get the record
  Record& record = S_->GetRecordW(key_, tag_next_, PK_type::getName());

  // Initialize the data
  if (!record.initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      Errors::Message message;
      message << PK_type::getName() << " has no initial condition parameter list.";
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
    markChangedSolutionPK(tag_next_);
  }

  PK_type::initialize();
};


// -----------------------------------------------------------------------------
// Called after all PKs have finished advancing a step successfully
// -----------------------------------------------------------------------------
template <class PK_type>
void
PK_Physical_Default<PK_type>::commitStep(double t_old, double t_new, const Tag& tag_next)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state @ " << tag_next << std::endl;

  // copy over primary variable
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(key_, tag_current, tag_next, *S_);

  PK_type::commitStep(t_old, t_new, tag_next);
}


// -----------------------------------------------------------------------------
// Called if any PKs have failed to advance a step successfully
// -----------------------------------------------------------------------------
template <class PK_type>
void
PK_Physical_Default<PK_type>::failStep(double t_old, double t_new, const Tag& tag_next)
{
  AMANZI_ASSERT(tag_next == tag_next_ || tag_next == Tags::NEXT);
  Tag tag_current = tag_next == tag_next_ ? tag_current_ : Tags::CURRENT;
  PKHelpers::assign(key_, tag_next, tag_current, *S_);

  PK_type::failStep(t_old, t_new, tag_next);
}


// -----------------------------------------------------------------------------
// Ensures the step size is smaller than max_valid_change
// -----------------------------------------------------------------------------
template <class PK_type>
bool
PK_Physical_Default<PK_type>::isValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Validating time step." << std::endl;

  if (max_valid_change_ > 0.0) {
    const CompositeVector& var_new = S_->template Get<CompositeVector>(key_, tag_next_);
    const CompositeVector& var_old = S_->template Get<CompositeVector>(key_, tag_current_);
    CompositeVector dvar(var_new.getMap());
    dvar.update(-1., var_old, 1., var_new, 0.);
    double change = dvar.normInf();
    if (change > max_valid_change_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max primary variable change=" << change
                   << " > limit=" << max_valid_change_ << std::endl;
      return false;
    }
  }
  return PK_type::isValidStep();
}


// -----------------------------------------------------------------------------
//  Marks as changed
// -----------------------------------------------------------------------------
template <class PK_type>
void
PK_Physical_Default<PK_type>::markChangedSolutionPK(const Tag& tag)
{
  PKHelpers::changedEvaluatorPrimary(key_, tag, *S_);
}

} // namespace Amanzi
