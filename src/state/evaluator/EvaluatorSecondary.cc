/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "EvaluatorSecondary.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
EvaluatorSecondary::EvaluatorSecondary(Teuchos::ParameterList& plist)
    : vo_(Keys::cleanPListName(plist.name()), plist), plist_(plist),
      computed_once_(false)
{
  // process the plist for names and tags of the things this evaluator
  // calculates
  if (plist_.isParameter("names")) {
    auto names = plist_.get<Teuchos::Array<std::string>>("names");
    if (plist_.isParameter("tags")) {
      auto tags = plist_.get<Teuchos::Array<std::string>>("tags");
      if (names.size() != tags.size()) {
        Errors::Message message;
        message << "EvaluatorSecondary: " << Keys::cleanPListName(plist.name())
                << " has names and tags lists of different sizes!";
        throw(message);
      }

      int i = 0;
      for (auto name : names) {
        my_keys_.emplace_back(std::make_pair(name, tags[i]));
        ++i;
      }
    } else {
      auto tag = plist_.get<std::string>("tag");
      for (auto name : names) {
        my_keys_.emplace_back(std::make_pair(name, tag));
      }
    }
  } else {
    auto name = Keys::cleanPListName(plist.name());
    if (plist_.isParameter("tags")) {
      auto tags = plist_.get<Teuchos::Array<std::string>>("tags");
      for (auto tag : tags) {
        my_keys_.emplace_back(std::make_pair(name, tag));
      }
    } else {
      auto tag = plist_.get<std::string>("tag");
      my_keys_.emplace_back(std::make_pair(name, tag));
    }
  }
  if (my_keys_.size() == 0) {
    Errors::Message message;
    message << "EvaluatorSecondary: " << plist.name()
            << " processed no key-tag pairs.";
    throw(message);
  }

  // process the plist for dependencies
  if (plist_.isParameter("dependencies")) {
    auto deps = plist_.get<Teuchos::Array<std::string>>("dependencies");

    if (plist_.isParameter("dependency tags")) {
      Teuchos::Array<std::string> tags =
        plist_.get<Teuchos::Array<std::string>>("dependency tags");
      if (deps.size() != tags.size()) {
        Errors::Message message;
        message << "EvaluatorSecondary: " << my_keys_[0].first
                << " has dependency and tag lists of different sizes!";
        throw(message);
      }

      int i = 0;
      for (auto dep : deps) {
        dependencies_.emplace_back(std::make_pair(dep, tags[i]));
        ++i;
      }
    } else if (plist_.get<bool>("dependency tags are my tag", true)) {
      auto my_tag = my_keys_[0].second;
      for (auto dep : deps) {
        dependencies_.emplace_back(std::make_pair(dep, my_tag));
      }
    } else {
      Errors::Message message;
      message << "EvalutorSecondary for " << my_keys_[0].first
              << " was not provided its dependencies' tags.";
      throw(message);
    }
  }
}

Evaluator&
EvaluatorSecondary::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorSecondary* other_p =
      dynamic_cast<const EvaluatorSecondary*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

EvaluatorSecondary&
EvaluatorSecondary::operator=(const EvaluatorSecondary& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_keys_ == other.my_keys_);
    requests_ = other.requests_;
    deriv_requests_ = other.deriv_requests_;
  }
  return *this;
}

// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool
EvaluatorSecondary::Update(State& S, const Key& request)
{
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "SecondaryVariable " << my_keys_[0].first << " requested by "
              << request << "..." << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = !computed_once_;
  for (auto& dep : dependencies_) {
    update |=
      S.GetEvaluator(dep.first, dep.second)
        .Update(S, Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second));
  }

  if (update) {
    computed_once_ = true;
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "... updating " << my_keys_[0].first << " value."
                << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    Update_(S);
    requests_.clear();
    requests_.insert(request);
    return true;
  } else {
    // Otherwise, see if we have filled this request already.
    if (requests_.find(request) == requests_.end()) {
      requests_.insert(request);
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << my_keys_[0].first
                  << "... has changed, but no need to update. " << std::endl;
      }
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << my_keys_[0].first << "... has not changed. " << std::endl;
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
bool
EvaluatorSecondary::UpdateDerivative(State& S, const Key& requestor,
                                     const Key& wrt_key, const Key& wrt_tag)
{
  if (!IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
    Errors::Message msg;
    msg << "EvaluatorSecondary (" << my_keys_[0].first << "," << my_keys_[0].second << ") is not differentiable with respect to (" << wrt_key << "," << wrt_tag << ").";
    throw(msg);
  }

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "SecondaryVariable Derivative d" << my_keys_[0].first << "_d" << wrt_key
              << " requested by " << requestor << "..." << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = false;

  // -- must update if our our dependencies have changed, as these affect the
  // partial derivatives
  Key my_request = Key{ "d" } +
                   Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second) +
                   "_d" + Keys::getKeyTag(wrt_key, wrt_tag);
  //  update |= Update(S, my_request);

  // -- must update if any of our dependencies or our dependencies' derivatives have changed
  for (auto& dep : dependencies_) {
    update |= S.GetEvaluator(dep.first, dep.second).Update(S, my_request);
    
    if (S.GetEvaluator(dep.first, dep.second)
          .IsDependency(S, wrt_key, wrt_tag)) {
      update |= S.GetEvaluator(dep.first, dep.second)
                  .UpdateDerivative(S, my_request, wrt_key, wrt_tag);
    }
  }

  // Do the update
  auto request = std::make_tuple(wrt_key, wrt_tag, requestor);
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "... updating." << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateDerivative_(S, wrt_key, wrt_tag);
    deriv_requests_.clear();
    deriv_requests_.insert(request);
    return true;
  } else {
    // Otherwise, simply service the request
    if (deriv_requests_.find(request) == deriv_requests_.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "... not updating but new to this request." << std::endl;
      }
      deriv_requests_.insert(request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "... has not changed." << std::endl;
      }
      return false;
    }
  }
}

inline bool
EvaluatorSecondary::IsDependency(const State& S, const Key& key,
                                 const Key& tag) const
{
  if (std::find(dependencies_.begin(),
                dependencies_.end(),
                std::make_pair(key, tag)) != dependencies_.end()) {
    return true;
  } else {
    for (auto& dep : dependencies_) {
      if (S.GetEvaluator(dep.first, dep.second).IsDependency(S, key, tag)) {
        return true;
      }
    }
  }
  return false;
}

inline bool
EvaluatorSecondary::ProvidesKey(const Key& key, const Key& tag) const
{
  return std::find(my_keys_.begin(),
                   my_keys_.end(),
                   std::make_pair(key, tag)) != my_keys_.end();
}

std::string
EvaluatorSecondary::WriteToString() const
{
  std::stringstream result;
  for (const auto& key : my_keys_) {
    result << key.first << ":" << (key.second.empty() ? "< >" : key.second) << ",";
  }
  result << std::endl << "  Type: secondary: " << name() << std::endl;
  for (const auto& dep : dependencies_) {
    result << "  Dep: " << dep.first << "," << dep.second << std::endl;
  }
  result << std::endl;
  return result.str();
}


void
EvaluatorSecondary::EnsureCompatibility_Flags_(State& S)
{
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  for (auto keytag : my_keys_) {
    S.GetRecordSet(keytag.first).set_io_vis(io_my_key);
    S.GetRecordSet(keytag.first).set_io_checkpoint(checkpoint_my_key);
  }
}


} // namespace Amanzi
