/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"

#define DEBUG_FLAG 0

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------
PrimaryVariableFieldEvaluator::PrimaryVariableFieldEvaluator(
    Teuchos::ParameterList& plist) :
    FieldEvaluator(plist) {
  type_ = EvaluatorType::PRIMARY;
  my_key_ = plist_.get<std::string>("evaluator name");
}

PrimaryVariableFieldEvaluator::PrimaryVariableFieldEvaluator(
        const PrimaryVariableFieldEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_) {
}

// ---------------------------------------------------------------------------
// Virtual copy constructor.
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator> PrimaryVariableFieldEvaluator::Clone() const {
  return Teuchos::rcp(new PrimaryVariableFieldEvaluator(*this));
}

void PrimaryVariableFieldEvaluator::operator=(const FieldEvaluator& other) {
#if ENABLE_DBC
  const PrimaryVariableFieldEvaluator* other_p =
    dynamic_cast<const PrimaryVariableFieldEvaluator*>(&other);
  AMANZI_ASSERT(other_p != NULL);
#else
  const PrimaryVariableFieldEvaluator* other_p =
    static_cast<const PrimaryVariableFieldEvaluator*>(&other);
#endif
  AMANZI_ASSERT(my_key_ == other_p->my_key_);

  requests_ = other_p->requests_;
  deriv_requests_ = other_p->deriv_requests_;
}


// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// by a field.
// ---------------------------------------------------------------------------
bool PrimaryVariableFieldEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Primary Variable " << my_key_ << " requested by " << request << std::endl;
  }

  if (requests_.find(request) == requests_.end()) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "  Has changed... " << std::endl;
    }
    requests_.insert(request);
    return true;

  } else {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "  Has not changed... " << std::endl;
    }
    return false;
  }
}


// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// ---------------------------------------------------------------------------
bool PrimaryVariableFieldEvaluator::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Primary Variable d" << my_key_ << "d" << wrt_key
          << " requested by " << request << std::endl;
  }

  std::pair<Key,Key> deriv_request(request, wrt_key);
  if (deriv_requests_.find(deriv_request) == deriv_requests_.end()) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "  Has changed..." << std::endl;
    }
    deriv_requests_.insert(deriv_request);
    return true;
  } else {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "  Has not changed..." << std::endl;
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
void PrimaryVariableFieldEvaluator::SetFieldAsChanged(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // clear cache
  requests_.clear();
  deriv_requests_.clear();

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Primary field \"" << vo_->color("gree") << my_key_ 
               << vo_->reset() << "\" marked as changed." << std::endl;
  }
}


bool PrimaryVariableFieldEvaluator::IsDependency(const Teuchos::Ptr<State>& S,
        Key key) const {
  return false;
}


bool PrimaryVariableFieldEvaluator::ProvidesKey(Key key) const {
  return (key == my_key_);
}

std::string
PrimaryVariableFieldEvaluator::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: primary" << std::endl
         << std::endl;
  return result.str();
}

} // namespace
