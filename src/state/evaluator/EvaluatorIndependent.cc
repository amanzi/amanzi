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
#include "EvaluatorIndependent.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependent::EvaluatorIndependent(Teuchos::ParameterList& plist) :
    time_(0.),
    computed_once_(false),
    my_key_(plist.name()),
    my_tag_(plist.get<std::string>("tag","")),
    temporally_variable_(!plist.get<bool>("constant in time", false)),
    plist_(plist),
    vo_(plist.name(), plist) {}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void EvaluatorIndependent::EnsureCompatibility(State& S) {
  // Require the field and claim ownership.
  S.Require<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_, my_key_);
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
bool EvaluatorIndependent::Update(State& S, const Key& request) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (!computed_once_) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Independent field \"" << my_key_ << "\" requested by " << request
            << " is updating for the first time." << std::endl;
    }

    // field DOES have to be computed at least once, even if it never changes.
    Update_(S);
    computed_once_ = true;
    return true;
  }

  if (temporally_variable_ && (S.time() != time_)) { // field is not current, update and clear requests
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Independent field \"" << my_key_ << "\" requested by " << request 
                 << " is updating." << std::endl;
    }
    Update_(S);
    requests_.clear();
    requests_.insert(request);
    return true;
  } else { // field is current, see if we have provided this request previously
    if (requests_.find(request) == requests_.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "Independent field \"" << vo_.color("green") << my_key_ 
                   << vo_.reset() << "\" requested by " << request << " has changed." << std::endl;
      }
      requests_.insert(request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "Independent field \"" << my_key_ << "\" requested by " << request 
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
bool EvaluatorIndependent::UpdateDerivative(State& S,
        const Key& request, const Key& wrt_key) {

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_.getOSTab();
    *vo_.os() << "INDEPENDENT Variable derivative requested by " << request
          << " has not changed." << std::endl;
  }
  return false;  // no derivatives, though this should never be called
}


bool EvaluatorIndependent::IsDependency(const State& S, const Key& key) const {
  return false;
}

bool EvaluatorIndependent::ProvidesKey(const Key& key) const {
  return key == my_key_;
}



// ---------------------------------------------------------------------------
// String representation of this evaluator
// ---------------------------------------------------------------------------
std::string
EvaluatorIndependent::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: independent" << std::endl
         << std::endl;
  return result.str();
}

} // namespace
