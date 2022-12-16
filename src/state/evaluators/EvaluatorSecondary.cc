/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Default field evaluator base class.  A Evaluator is a node in the
  dependency graph.
*/

#include "EvaluatorSecondary.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
EvaluatorSecondary::EvaluatorSecondary(Teuchos::ParameterList& plist)
  : nonlocal_dependencies_(false),
    updated_once_(false),
    plist_(plist),
    vo_(Keys::cleanPListName(plist.name()), plist)
{
  type_ = EvaluatorType::SECONDARY;

  // process the plist for names and tags of the things this evaluator calculates
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
        my_keys_.emplace_back(KeyTag{ name, Tag{ tags[i] } });
        ++i;
      }
    } else {
      auto tag = make_tag(plist_.get<std::string>("tag"));
      for (auto name : names) { my_keys_.emplace_back(KeyTag{ name, Tag{ tag } }); }
    }
  } else {
    auto name = Keys::cleanPListName(plist.name());
    if (plist_.isParameter("tags")) {
      auto tags = plist_.get<Teuchos::Array<std::string>>("tags");
      for (auto tag : tags) { my_keys_.emplace_back(KeyTag{ name, Tag{ tag } }); }
    } else {
      auto tag = make_tag(plist_.get<std::string>("tag"));
      my_keys_.emplace_back(KeyTag{ name, Tag{ tag } });
    }
  }
  if (my_keys_.size() == 0) {
    Errors::Message msg;
    msg << "EvaluatorSecondary: " << plist.name() << " processed no key-tag pairs.";
    throw(msg);
  }

  // process the plist for dependencies
  if (plist_.isParameter("dependencies")) {
    auto deps = plist_.get<Teuchos::Array<std::string>>("dependencies");

    if (plist_.isParameter("dependency tags")) {
      auto tags = plist_.get<Teuchos::Array<std::string>>("dependency tags");
      if (deps.size() != tags.size()) {
        Errors::Message message;
        message << "EvaluatorSecondary: " << my_keys_[0].first
                << " has dependency and tag lists of different sizes!";
        throw(message);
      }

      int i = 0;
      for (auto dep : deps) {
        dependencies_.insert(KeyTag(dep, make_tag(tags[i])));
        ++i;
      }
    } else {
      auto my_tag = my_keys_.front().second;
      auto dep_tag = Keys::readTag(plist_, my_tag);
      for (const auto& dep : deps) { dependencies_.insert(KeyTag(dep, dep_tag)); }
    }

  } else if (plist_.isParameter("dependency suffixes")) {
    auto dep_suff = plist_.get<Teuchos::Array<std::string>>("dependency suffixes");
    auto domain = Keys::getDomain(my_keys_.front().first);

    if (plist_.isParameter("dependency tags")) {
      auto tags = plist_.get<Teuchos::Array<std::string>>("dependency tags");
      if (dep_suff.size() != tags.size()) {
        Errors::Message message;
        message << "EvaluatorSecondary: " << my_keys_[0].first
                << " has dependency and tag lists of different sizes!";
        throw(message);
      }

      int i = 0;
      for (const auto& suffix : dep_suff) {
        dependencies_.insert(KeyTag(Keys::getKey(domain, suffix), make_tag(tags[i])));
        ++i;
      }
    } else {
      const auto& my_tag = my_keys_[0].second;
      auto dep_tag = Keys::readTag(plist_, my_tag);
      for (const auto& suffix : dep_suff) {
        dependencies_.insert(KeyTag(Keys::getKey(domain, suffix), dep_tag));
      }
    }
  }

  nonlocal_dependencies_ = plist_.get<bool>("includes non-rank-local dependencies", false);
}


Evaluator&
EvaluatorSecondary::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorSecondary* other_p = dynamic_cast<const EvaluatorSecondary*>(&other);
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


// ---------------------------------------------------------------------------
// Step 1 of graph checking -- Requires evaluators for the full dependency
// graph.
// ---------------------------------------------------------------------------
void
EvaluatorSecondary::EnsureEvaluators(State& S)
{
  // first make sure that we have evaluators for all of _our_ keys.  Note that
  // we may not, as not all of my_keys_ may be required.
  for (const auto& my_key_tag : my_keys_) {
    S.RequireEvaluator(my_key_tag.first, my_key_tag.second);
  }

  // next make sure that all dependencies have evaluators, and call
  // EnsureEvaluators recursively to fill out the dependency graph
  for (const auto& dep : dependencies_) {
    auto& dep_eval = S.RequireEvaluator(dep.first, dep.second);
    dep_eval.EnsureEvaluators(S);
  }
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
    *vo_.os() << "SecondaryVariable " << my_keys_[0].first << " requested by " << request
              << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = false;
  for (auto& dep : dependencies_) {
    update |= S.GetEvaluator(dep.first, dep.second)
                .Update(S, Keys::getKey(my_keys_[0].first, my_keys_[0].second));
  }

  if (!updated_once_) {
    // force always updating once
    update = true;
    updated_once_ = true;
  }

  // check if nonlocal for changes in offprocess dependencies
  if (nonlocal_dependencies_) {
    auto comm = get_comm_(S);
    if (comm != Teuchos::null) {
      int update_l = update;
      int update_g = 0;
      comm->MaxAll(&update_l, &update_g, 1);
      update |= update_g;
    }
  }

  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Updating " << my_keys_[0].first << " value... " << std::endl;
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
        *vo_.os() << my_keys_[0].first << " has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << my_keys_[0].first << " has not changed... " << std::endl;
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
EvaluatorSecondary::UpdateDerivative(State& S,
                                     const Key& requestor,
                                     const Key& wrt_key,
                                     const Tag& wrt_tag)
{
  AMANZI_ASSERT(IsDependency(S, wrt_key, wrt_tag));

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Algebraic Variable d" << my_keys_[0].first << "_d" << wrt_key << " requested by "
              << requestor << std::endl;
  }

  // If wrt_key is not a dependency, no need to differentiate.
  if (!IsDependency(S, wrt_key, wrt_tag)) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << wrt_key << " is not a dependency... " << std::endl;
    }
    return false;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = false;

  // -- must update if our our dependencies have changed, as these affect the
  // partial derivatives
  Key my_request = Keys::getDerivKey(my_keys_[0], KeyTag(wrt_key, wrt_tag));
  update |= Update(S, my_request);

  // -- must update if any of our dependencies' derivatives have changed
  for (auto& dep : dependencies_) {
    if (S.GetEvaluator(dep.first, dep.second).IsDependency(S, wrt_key, wrt_tag)) {
      update |=
        S.GetEvaluator(dep.first, dep.second).UpdateDerivative(S, my_request, wrt_key, wrt_tag);
    }
  }

  // check if nonlocal for changes in offprocess dependencies
  if (nonlocal_dependencies_) {
    auto comm = get_comm_(S);
    if (comm != Teuchos::null) {
      int update_l = update;
      int update_g = 0;
      comm->MaxAll(&update_l, &update_g, 1);
      update |= update_g;
    }
  }

  // Do the update
  DerivativeTriple request = std::make_tuple(wrt_key, wrt_tag, requestor);
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) { *vo_.os() << "  ... updating." << std::endl; }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateDerivative_(S, wrt_key, wrt_tag);
    deriv_requests_.clear();
    deriv_requests_.insert(request);
    return true;

  } else {
    // Otherwise, simply service the request
    if (deriv_requests_.find(request) == deriv_requests_.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "  ... not updating but new to this request." << std::endl;
      }
      deriv_requests_.insert(request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) { *vo_.os() << "  ... has not changed." << std::endl; }
      return false;
    }
  }
}


bool
EvaluatorSecondary::IsDependency(const State& S, const Key& key, const Tag& tag) const
{
  if (std::find(dependencies_.begin(), dependencies_.end(), std::make_pair(key, tag)) !=
      dependencies_.end()) {
    return true;
  } else {
    for (auto& dep : dependencies_) {
      if (S.GetEvaluator(dep.first, dep.second).IsDependency(S, key, tag)) { return true; }
    }
  }
  return false;
}


bool
EvaluatorSecondary::ProvidesKey(const Key& key, const Tag& tag) const
{
  std::pair<Key, Tag> my_key = std::make_pair(key, tag);
  return std::find(my_keys_.begin(), my_keys_.end(), my_key) != my_keys_.end();
}


std::string
EvaluatorSecondary::WriteToString() const
{
  std::stringstream result;
  for (const auto& key : my_keys_) { result << Keys::getKey(key.first, key.second); }
  result << std::endl << "  type: secondary" << std::endl;
  for (const auto& dep : dependencies_) {
    result << "  dep: " << dep.first << ", tag=\"" << dep.second.get() << "\"" << std::endl;
  }
  return result.str();
}


void
EvaluatorSecondary::EnsureCompatibility_Flags_(State& S)
{
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  for (auto keytag : my_keys_) {
    S.GetRecordW(keytag.first, keytag.second, keytag.first).set_io_vis(io_my_key);
    S.GetRecordW(keytag.first, keytag.second, keytag.first).set_io_checkpoint(checkpoint_my_key);
  }
}

} // namespace Amanzi
