/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#include "EvaluatorIndependent.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependent_::EvaluatorIndependent_(Teuchos::ParameterList& plist)
    : my_key_(Keys::cleanPListName(plist.name())),
      my_tag_(Keys::readTag(plist, "tag")),
      time_(0.0),
      temporally_variable_(!plist.get<bool>("constant in time", false)),
      computed_once_(false),
      plist_(plist),
      vo_(Keys::cleanPListName(plist.name()), plist) {
  type_ = EvaluatorType::INDEPENDENT;
}


// ---------------------------------------------------------------------------
// Assignment operator.
// ---------------------------------------------------------------------------
EvaluatorIndependent_& EvaluatorIndependent_::operator=(const EvaluatorIndependent_& other) {
  if (this !=& other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    time_ = other.time_;
    computed_once_ = other.computed_once_;
    temporally_variable_ = other.temporally_variable_;
    requests_ = other.requests_;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Virtual assignment operator.
// ---------------------------------------------------------------------------
Evaluator& EvaluatorIndependent_::operator=(const Evaluator& other) {
  if (this !=& other) {
    const EvaluatorIndependent_* other_p = dynamic_cast<const EvaluatorIndependent_*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void EvaluatorIndependent_::EnsureCompatibility(State& S) {
  S.GetRecordW(my_key_, my_tag_, my_key_).set_initialized();

  // check plist for vis or checkpointing control
  auto vis_check = std::string{"visualize " + my_key_};
  bool io_my_key = plist_.get<bool>(vis_check, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);

  auto chkp_check = std::string{"checkpoint " + my_key_};
  bool checkpoint_my_key = plist_.get<bool>(chkp_check, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);
}


// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// ---------------------------------------------------------------------------
bool EvaluatorIndependent_::Update(State& S, const Key& request) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (!computed_once_) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Independent field \"" << my_key_ << "\" requested by "
                << request << " is updating for the first time." << std::endl;
    }

    // field DOES have to be computed at least once, even if it never changes.
    Update_(S);
    requests_.insert(request);
    computed_once_ = true;
    if (temporally_variable_)
      time_ = S.get_time(my_tag_);
    return true;
  }

  if (temporally_variable_ && (S.get_time(my_tag_) != time_)) {
    // field is not current, update and clear requests
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Independent field \"" << my_key_ << "\" requested by "
                << request << " is updating." << std::endl;
    }
    Update_(S);
    if (temporally_variable_)
      time_ = S.get_time(my_tag_);
    requests_.clear();
    requests_.insert(request);
    return true;
  } else {
    // field is current, see if we have provided this request previously
    if (requests_.find(request) == requests_.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "Independent field \"" << vo_.color("green") << my_key_
                  << vo_.reset() << "\" requested by " << request
                  << " has changed." << std::endl;
      }
      requests_.insert(request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "Independent field \"" << my_key_ << "\" requested by "
                  << request << " has not changed." << std::endl;
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
bool EvaluatorIndependent_::UpdateDerivative(State& S, const Key& request,
                                             const Key& wrt_key,
                                             const Tag& wrt_tag) {

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_.getOSTab();
    *vo_.os() << "INDEPENDENT Variable derivative requested by " << request
              << " has not changed." << std::endl;
  }
  return false; // no derivatives, though this should never be called
}


bool EvaluatorIndependent_::IsDependency(const State& S, const Key& key, const Tag& tag) const {
  return false;
}


bool EvaluatorIndependent_::ProvidesKey(const Key& key, const Tag& tag) const {
  return (key == my_key_) && (tag == my_tag_);
}


bool EvaluatorIndependent_::IsDifferentiableWRT(const State& S, const Key& wrt_key,
                                                const Tag& wrt_tag) const {
  return false;
}


// void EvaluatorIndependent_::EnsureCompatibleDerivative(
//     State& S, const Key& wrt_key, const Tag& wrt_tag) {
//   Errors::Message msg("Independent Variables are not differentiable");
//   throw(msg);
// }


// ---------------------------------------------------------------------------
// String representation of this evaluator
// ---------------------------------------------------------------------------
std::string EvaluatorIndependent_::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: independent" << std::endl;
  return result.str();
}



} // namespace Amanzi
