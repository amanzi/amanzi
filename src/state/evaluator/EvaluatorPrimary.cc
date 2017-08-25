/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Arcos

License: BSD
Author: Ethan Coon

An evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#include "UniqueHelpers.hh"
#include "EvaluatorPrimary.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------
EvaluatorPrimary::EvaluatorPrimary(Teuchos::ParameterList& plist) :
    my_key_(plist.name()),
    my_tag_(plist.get<std::string>("tag", "")),
    vo_(plist.name(), plist) {}

// ---------------------------------------------------------------------------
// Virtual copy constructor.
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator> EvaluatorPrimary::Clone() const {
  return Teuchos::rcp(new EvaluatorPrimary(*this));
}

// ---------------------------------------------------------------------------
// Virtual assignment operator.
// ---------------------------------------------------------------------------
Evaluator& EvaluatorPrimary::operator=(const Evaluator& other) {
  if (this != &other) {
    const EvaluatorPrimary* other_p =
        dynamic_cast<const EvaluatorPrimary*>(&other);
    ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

// ---------------------------------------------------------------------------
// Assignment operator.
// ---------------------------------------------------------------------------
EvaluatorPrimary& EvaluatorPrimary::operator=(const EvaluatorPrimary& other) {
  if (this != &other) {
    ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
    deriv_requests_ = other.deriv_requests_;
  }
  return *this;
}

// ---------------------------------------------------------------------------
// Lazy evaluation of the evaluator.
//
// Updates the data, if needed.  Returns true if the value of the data has
// changed since the last request for an update.
// ---------------------------------------------------------------------------
bool EvaluatorPrimary::Update(State& S,
        const Key& request) {
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
bool EvaluatorPrimary::UpdateDerivative(State& S,
        const Key& request, const Key& wrt_key) {
  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Primary Variable d" << my_key_ << "d" << wrt_key
          << " requested by " << request << std::endl;
  }

  std::pair<Key,Key> deriv_request(request, wrt_key);
  if (deriv_requests_.find(deriv_request) == deriv_requests_.end()) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  Has changed..." << std::endl;
    }
    deriv_requests_.insert(deriv_request);
    return true;
  } else {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "  Has not changed..." << std::endl;
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
void EvaluatorPrimary::SetChanged() {
  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Primary field \"" << vo_.color("gree") << my_key_ 
               << vo_.reset() << "\" marked as changed." << std::endl;
  }

  // clear cache
  requests_.clear();
  deriv_requests_.clear();
}


// Nothing is a dependency of a primary variable.
bool EvaluatorPrimary::IsDependency(const State& S,
        const Key& key) const {
  return false;
}


bool EvaluatorPrimary::ProvidesKey(const Key& key) const {
  return (key == my_key_);
}

std::string
EvaluatorPrimary::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: primary" << std::endl
         << std::endl;
  return result.str();
}

} // namespace
