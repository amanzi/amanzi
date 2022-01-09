/* -------------------------------------------------------------------------
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  An evaluator with no dependencies solved for by a PK.
*/

#include "EvaluatorPrimary.hh"
#include "StateDefs.hh"
#include "UniqueHelpers.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------
EvaluatorPrimary_::EvaluatorPrimary_(Teuchos::ParameterList& plist)
    : my_key_(Keys::cleanPListName(plist.name())),
      my_tag_(make_tag(plist.get<std::string>("tag", ""))),
      vo_(Keys::cleanPListName(plist.name()), plist)
{
  type_ = EvaluatorType::PRIMARY;
}


// ---------------------------------------------------------------------------
// Assignment operator.
// ---------------------------------------------------------------------------
EvaluatorPrimary_& EvaluatorPrimary_::operator=(const EvaluatorPrimary_& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
    deriv_requests_ = other.deriv_requests_;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Virtual assignment operator.
// ---------------------------------------------------------------------------
Evaluator& EvaluatorPrimary_::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorPrimary_ *other_p = dynamic_cast<const EvaluatorPrimary_*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Lazy evaluation of the evaluator.
//
// Updates the data, if needed.  Returns true if the value of the data has
// changed since the last request for an update.
// ---------------------------------------------------------------------------
bool EvaluatorPrimary_::Update(State& S, const Key& request)
{
  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Primary Variable " << my_key_ << " requested by " << request << std::endl;
  }

  if (requests_.find(request) == requests_.end()) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  Has changed... " << std::endl;
    }
    requests_.insert(request);
    return true;

  } else {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  Has not changed... " << std::endl;
    }
    return false;
  }
}


// ---------------------------------------------------------------------------
// Lazy evaluation of derivatives of evaluator.
//
// Updates the derivative, if needed.  Returns true if the value of the
// derivative with respect to wrt_key has changed since the last request for
// an update.
// ---------------------------------------------------------------------------
bool EvaluatorPrimary_::UpdateDerivative(State& S, const Key& request,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag)
{
  // enforce the contract: all calls to this must either have wrt_key,wrt_tag
  // as the key provided (for primary evaluators) or as a dependency (for
  // secondary evaluators)
  AMANZI_ASSERT(ProvidesKey(wrt_key, wrt_tag));

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Primary Variable " << my_key_ << ":" << my_tag_.get()
              << " derivative with respect to " << wrt_key << ":" << wrt_tag.get()
              << " requested by " << request;
  }

  if (deriv_requests_.size() == 0) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  ... updating." << std::endl;
    }
    // derivative with respect to me, lazy calculate, answer is 1
    UpdateDerivative_(S);
  }

  if (deriv_requests_.find(std::make_tuple(wrt_key, wrt_tag, request)) == deriv_requests_.end()) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  ... not updating but new to this request." << std::endl;
    }
    deriv_requests_.insert(std::make_tuple(wrt_key, wrt_tag, request));
    return true;
  } else {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  ... has not changed." << std::endl;
    }
    return false;
  }
}


// ---------------------------------------------------------------------------
// How a PK informs this leaf of the tree that it has changed.
//
// Effectively this simply tosses the request history, so that the next
// requests will say this has changed.
// ---------------------------------------------------------------------------
void EvaluatorPrimary_::SetChanged()
{
  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Primary field \"" << vo_.color("gree") << my_key_
              << vo_.reset() << "\" marked as changed." << std::endl;
  }

  // clear cache
  requests_.clear();

  // note the derivative is always 1, no need to re-calculate.
}


// Nothing is a dependency of a primary variable.
bool EvaluatorPrimary_::IsDependency(const State& S, const Key& key, const Tag& tag) const {
  return false;
}


bool EvaluatorPrimary_::ProvidesKey(const Key& key, const Tag& tag) const {
  return (key == my_key_) && (tag == my_tag_);
}


std::string EvaluatorPrimary_::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl << "  Type: primary" << std::endl;
  return result.str();
}

} // namespace Amanzi
