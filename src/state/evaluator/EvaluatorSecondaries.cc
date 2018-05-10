/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default field evaluator base class.  A Evaluator is a node in the dependency
graph.

------------------------------------------------------------------------- */
#include "EvaluatorSecondaries.hh"

#include "EvaluatorPrimary.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
EvaluatorSecondaries::EvaluatorSecondaries(Teuchos::ParameterList &plist)
    : vo_(Keys::cleanPListName(plist.name()), plist), plist_(plist) {
  // process the plist
  if (plist_.isParameter("names")) {
    auto names = plist_.get<Teuchos::Array<std::string>>("names");
    if (plist_.isParameter("tags")) {
      auto tags = plist_.get<Teuchos::Array<std::string>>("tags");
      if (names.size() != tags.size()) {
        Errors::Message message;
        message << "EvaluatorSecondaries: "
                << Keys::cleanPListName(plist.name())
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
  }

  if (plist_.isParameter("dependencies")) {
    Teuchos::Array<std::string> deps =
        plist_.get<Teuchos::Array<std::string>>("dependencies");
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
    } else if (plist_.get<bool>("dependency tags are my tag", false)) {
      auto my_tag = plist_.get<std::string>("evaluator tag");
      for (auto dep : deps) {
        dependencies_.emplace_back(std::make_pair(dep, my_tag));
      }
    } else {
      Errors::Message message;
      message << "EvalutorSecondaries for " << my_keys_[0].first
              << " was not provided its dependencies' tags.";
      throw(message);
    }
  }
}

Evaluator &EvaluatorSecondaries::operator=(const Evaluator &other) {
  if (this != &other) {
    const EvaluatorSecondaries *other_p =
        dynamic_cast<const EvaluatorSecondaries *>(&other);
    ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

EvaluatorSecondaries &EvaluatorSecondaries::
operator=(const EvaluatorSecondaries &other) {
  if (this != &other) {
    ASSERT(my_keys_ == other.my_keys_);
    requests_ = other.requests_;
    deriv_requests_ = other.deriv_requests_;
  }
  return *this;
}

// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool EvaluatorSecondaries::Update(State &S, const Key &request) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "SecondariesVariable " << my_keys_[0].first << " requested by "
              << request << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our
  // dependencies.
  bool update = false;
  for (auto &dep : dependencies_) {
    update |= S.GetEvaluator(dep.first, dep.second)
                  .Update(S, Keys::getKeyTag(my_keys_[0].first,
                                               my_keys_[0].second));
  }

  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Updating " << my_keys_[0].first << " value... "
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
                  << " has changed, but no need to update... " << std::endl;
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
bool EvaluatorSecondaries::UpdateDerivative(State &S, const Key &requestor,
                                            const Key &wrt_key,
                                            const Key &wrt_tag) {
  ASSERT(IsDependency(S, wrt_key, wrt_tag));

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Algebraic Variable d" << my_keys_[0].first << "_d" << wrt_key
              << " requested by " << requestor << std::endl;
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
  Key my_request = Key{"d"} +
                   Keys::getKeyTag(my_keys_[0].first, my_keys_[0].second) +
                   "_d" + Keys::getKeyTag(wrt_key, wrt_tag);
  update |= Update(S, my_request);

  // -- must update if any of our dependencies' derivatives have changed
  for (auto &dep : dependencies_) {
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
      *vo_.os() << "  ... updating." << std::endl;
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
        *vo_.os() << "  ... not updating but new to this request." << std::endl;
      }
      deriv_requests_.insert(request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "  ... has not changed." << std::endl;
      }
      return false;
    }
  }
}

// ---------------------------------------------------------------------------
// Updates the field value in state S.
// ---------------------------------------------------------------------------
void EvaluatorSecondaries::Update_(State& S) {
  // pull my variables out of state
  std::vector<Teuchos::Ptr<CompositeVector> > myfields;
  for (auto& key : my_keys_) {
    myfields.push_back(S.GetPtrW<CompositeVector>(key.first, key.second, key.first).ptr());
  }

  // call the evaluate method
  Evaluate_(S, myfields);
}


// ---------------------------------------------------------------------------
// Updates the total derivative d(my_field)/d(wrt_field) in state S.
//
// (using p to indicate partial derivative...)
// d(my_field)/d(wrt_field) = p(my_field)/p(wrt_field) + sum_(dep) d(dep)/d(wrt_field)
// ---------------------------------------------------------------------------
void EvaluatorSecondaries::UpdateDerivative_(State& S,
                                             const Key& wrt_key,
                                             const Key& wrt_tag) {
  
  std::vector<Teuchos::Ptr<CompositeVector> > dmys;
  for (auto& my_key : my_keys_) {

    // derivative of my_key.first(my_key.second) with respect to wrt_key(wrt_tag)
    
    Key dmy_key = std::string("d")+my_key.first+std::string("_d")+wrt_key;

    if (S.HasData(dmy_key, my_key.second )) {
      // Get the field...
      Teuchos::RCP<CompositeVector> dmy = S.GetPtrW<CompositeVector>(dmy_key, my_key.second, my_key.first);
      dmy->PutScalar(0.0);
      dmys.push_back(dmy.ptr());
    } else {
      // or create the field.  Note we have to do extra work that is normally
      // done by State in initialize.
      CompositeVectorSpace& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key.first, my_key.second);
      CompositeVectorSpace& new_fac = S.Require<CompositeVector,CompositeVectorSpace>(dmy_key, my_key.second, my_key.first);
      new_fac.Update(my_fac);
      Teuchos::RCP<CompositeVector> dmy = Teuchos::rcp(new CompositeVector(new_fac));
      S.SetPtr<CompositeVector>(dmy_key, my_key.second, my_key.first, dmy);
      S.GetRecordW(dmy_key, my_key.first).set_initialized();
      S.GetRecordW(dmy_key, my_key.first).set_io_vis(plist_.get<bool>("visualize derivative", false));
      S.GetRecordW(dmy_key, my_key.first).set_io_checkpoint(plist_.get<bool>("checkpoint derivative", false));

      dmy->PutScalar(0.0);
      dmys.push_back(dmy.ptr());
    }
  }

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto& dep : dependencies_) {

    // -- allocate a tmp set of vectors for dF_i/ddep
    std::vector<Teuchos::Ptr<CompositeVector> > tmp(dmys.size());
    for (int i=0; i!=dmys.size(); ++i) {
      tmp[i] = Teuchos::ptr(new CompositeVector(*dmys[i]));
    }

    if ((wrt_key == dep.first)&&(wrt_tag==dep.second)) {
      
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      for (int i=0; i!=dmys.size(); ++i) {
        // partial F_i / partial x
        dmys[i]->Update(1.0, *tmp[i], 1.0);
      }
      
    } else if (S.GetEvaluator(dep.first, dep.second)->IsDependency(S, wrt_key, wrt_tag)) {
      //tmp = partial F_i / partial dep
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);
      
      for (int i=0; i!=dmys.size(); ++i) {
        // -- ddep/dx
        Key ddep_key = std::string("d")+dep.first+std::string("_d")+wrt_key;
        const CompositeVector& ddep = S.Get<CompositeVector>(ddep_key, dep.second);
        // partial F_i / partial dep * ddep/dx
        dmys[i]->Multiply(1.0, ddep, *tmp[i], 1.0);
      }
    }

    // clean up tmp as it goes out of scope
    for (int i=0; i!=dmys.size(); ++i) {
      delete tmp[i].get();
    }
  }
  // // if (check_derivative_) {
  // //   CheckDerivative_(S, wrt_key, wrt_tag);
  // // }
  
}

  

inline bool EvaluatorSecondaries::IsDependency(const State &S, const Key &key,
                                               const Key &tag) const {
  if (std::find(dependencies_.begin(), dependencies_.end(),
                std::make_pair(key, tag)) != dependencies_.end()) {
    return true;
  } else {
    for (auto &dep : dependencies_) {
      if (S.GetEvaluator(dep.first, dep.second).IsDependency(S, key, tag)) {
        return true;
      }
    }
  }
  return false;
}


bool EvaluatorSecondaries::ProvidesKey(const Key& key, const Key& tag) const {
  return std::find(my_keys_.begin(), my_keys_.end(), std::make_pair(key,tag)) != my_keys_.end();
}

void EvaluatorSecondaries::EnsureCompatibility(State& S) {

  Teuchos::RCP<CompositeVectorSpace> master_fac;

  for (auto& key : my_keys_) {
    // Ensure my field exists, and claim ownership.
    CompositeVectorSpace& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(key.first, key.second, key.first);

    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>("visualize", true);
    S.GetRecordW(key.first, key.second).set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
    S.GetRecordW(key.first, key.second).set_io_checkpoint(checkpoint_my_key);

    // Select a "master factory" to ensure commonality of all of my factories.
    if (my_fac.Mesh() != Teuchos::null) {
      if (master_fac == Teuchos::null) {
        // no master yet, take this one
        master_fac = Teuchos::rcpFromRef(my_fac);
      } else {
        // existing master, ensure it includes these requirements
        master_fac->Update(my_fac);
      }
    }
  }

  if (master_fac == Teuchos::null) {
    // No requirements have been set, so we'll have to hope they get set by an
    // evaluator that depends upon this evaluator.
  } else {
    // Create an unowned factory to check my dependencies.
    // Teuchos::RCP<CompositeVectorSpace> dep_fac =
    //     Teuchos::rcp(new CompositeVectorSpace(*master_fac));
    // dep_fac->SetOwned(false);
    
    CompositeVectorSpace dep_fac(*master_fac);
    dep_fac.SetOwned(false);
    
    // Loop over my keys, ensuring they meet the requirements of the master.
    for (auto& key : my_keys_ ) {
      CompositeVectorSpace& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(key.first, key.second);
      my_fac.Update(dep_fac);
    }

    // Loop over my dependencies, ensuring they meet my requirements.
    for (auto& key : dependencies_) {
      CompositeVectorSpace& fac = S.Require<CompositeVector,CompositeVectorSpace>(key.first, key.second);
      fac.Update(dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (auto& key : dependencies_) {
      S.RequireEvaluator(key.first, key.second)->EnsureCompatibility(S);
    }
    
  }
}
  

std::string EvaluatorSecondaries::WriteToString() const {
  std::stringstream result;
  for (const auto &key : my_keys_) {
    result << key.first << ":" << key.second << ",";
  }
  result << std::endl << "  Type: secondary" << std::endl;
  for (const auto &dep : dependencies_) {
    result << "  Dep: " << dep.first << "," << dep.second << std::endl;
  }
  result << std::endl;
  return result.str();
}

} // namespace Amanzi
