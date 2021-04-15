/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#include "CompositeVectorFunctionFactory.hh"
#include "independent_variable_field_evaluator.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
IndependentVariableFieldEvaluator::IndependentVariableFieldEvaluator(Teuchos::ParameterList& plist) :
    FieldEvaluator(plist),
    time_(0.),
    temporally_variable_(true),
    computed_once_(false) {

  type_ = EvaluatorType::INDEPENDENT;
  my_key_ = plist_.get<std::string>("evaluator name");
  temporally_variable_ = !plist_.get<bool>("constant in time", false);
}


// ---------------------------------------------------------------------------
// Copy constructor
// ---------------------------------------------------------------------------
IndependentVariableFieldEvaluator::IndependentVariableFieldEvaluator(const IndependentVariableFieldEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_),
    temporally_variable_(other.temporally_variable_),
    time_(other.time_),
    computed_once_(other.computed_once_) {}


// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
void IndependentVariableFieldEvaluator::operator=(const FieldEvaluator& other) {
#if ENABLE_DBC
  const IndependentVariableFieldEvaluator* other_p =
      dynamic_cast<const IndependentVariableFieldEvaluator*>(&other);
  AMANZI_ASSERT(other_p != NULL);
#else
  const IndependentVariableFieldEvaluator* other_p =
      static_cast<const IndependentVariableFieldEvaluator*>(&other);
#endif
  AMANZI_ASSERT(my_key_ == other_p->my_key_);

  time_ = other_p->time_;
  computed_once_ = other_p->computed_once_;
  temporally_variable_ = other_p->temporally_variable_;
  requests_ = other_p->requests_;
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void IndependentVariableFieldEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Require the field and claim ownership.
  const auto& fac = S->RequireField(my_key_, my_key_);

  // set not owned -- this allows subsequent calls to add to my components
  fac->SetOwned(false);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
}

// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// ---------------------------------------------------------------------------
bool IndependentVariableFieldEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (!computed_once_) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Independent field \"" << my_key_ << "\" requested by " << request
            << " is updating for the first time." << std::endl;
    }

    // field DOES have to be computed at least once, even if it never changes.
    time_ = S->time();
    UpdateField_(S);
    computed_once_ = true;
    requests_.insert(request);
    return true;
  }

  if (temporally_variable_ && (S->time() != time_)) { // field is not current, update and clear requests
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Independent field \"" << my_key_ << "\" evaluated previously at time = "
                 << time_ << " requested by " << request 
                 << " is updating at time = " << S->time() << std::endl;
    }
    time_ = S->time();
    UpdateField_(S);
    requests_.clear();
    requests_.insert(request);
    return true;
  } else { // field is current, see if we have provided this request previously
    if (requests_.find(request) == requests_.end()) {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "Independent field \"" << vo_->color("green") << my_key_ 
                   << vo_->reset() << "\" requested by " << request << " has changed." << std::endl;
      }
      requests_.insert(request);
      return true;
    } else {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "Independent field \"" << my_key_ << "\" requested by " << request 
                   << " has not changed." << std::endl;
      }
      return false;
    }
  }
}


// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
bool IndependentVariableFieldEvaluator::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "INDEPENDENT Variable derivative requested by " << request
          << " has not changed." << std::endl;
  }
  return false;  // no derivatives, though this should never be called
}


bool IndependentVariableFieldEvaluator::IsDependency(const Teuchos::Ptr<State>& S, Key key) const {
  return false;
}

bool IndependentVariableFieldEvaluator::ProvidesKey(Key key) const { return key == my_key_; }


// ---------------------------------------------------------------------------
// String representation of this evaluator
// ---------------------------------------------------------------------------
std::string
IndependentVariableFieldEvaluator::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: independent" << std::endl
         << std::endl;
  return result.str();
}

} // namespace
