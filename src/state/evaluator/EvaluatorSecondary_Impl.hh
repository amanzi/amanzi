/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default field evaluator base class.  A Evaluator is a node in the dependency
graph.

------------------------------------------------------------------------- */


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
EvaluatorSecondary<Data_t, DataFactory_t>::EvaluatorSecondary(Teuchos::ParameterList& plist) :
    my_key_(plist.name()),
    my_tag_(plist.get<std::string>("tag","")),
    vo_(plist.name(), plist),
    plist_(plist),
    check_derivative_(plist.get<bool>("check derivatives", false))
{
  if (plist_.isParameter("dependencies")) {
    Teuchos::Array<std::string> deps =
        plist_.get<Teuchos::Array<std::string> >("dependencies");
    if (plist_.isParameter("dependency tags")) {
      Teuchos::Array<std::string> tags =
          plist_.get<Teuchos::Array<std::string> >("dependency tags");
      if (deps.size() != tags.size()) {
        Errors::Message message;
        message << "EvaluatorSecondary: " << my_key_ << " has dependency and tag lists of different sizes!";
        throw(message);
      }

      int i=0;
      for (auto dep : deps) {
        dependencies_.emplace_back(std::make_pair(dep, tags[i]));
        ++i;
      }
    } else if (plist_.get<bool>("dependency tags are my tag", false)) {
      for (auto dep : deps) {
        dependencies_.emplace_back(std::make_pair(dep, my_tag_));
      }
    } else {
      Errors::Message message;
      message << "EvalutorSecondary for " << my_key_ << " was not provided its dependencies' tags.";
      throw(message);
    }
  }
}


// -----------------------------------------------------------------------------
// Assignment operators
// -----------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
Evaluator&
EvaluatorSecondary<Data_t, DataFactory_t>::operator=(const Evaluator& other) {
  if (this != &other) {
    const EvaluatorSecondary* other_p =
        dynamic_cast<const EvaluatorSecondary*>(&other);
    ASSERT(other_p != NULL);
    ASSERT(my_key_ == other_p->my_key_);

    my_tag_ = other_p->my_tag_;
    requests_ = other_p->requests_;
    deriv_requests_ = other_p->deriv_requests_;
    dependencies_ = other_p->dependencies_;
    check_derivative_ = other_p->check_derivative_;
  }
  return *this;
}

// -----------------------------------------------------------------------------
// Answers the question, has this data changed since it was last requested
// by a given requestor.  Updates the data if needed.
// -----------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
bool EvaluatorSecondary<Data_t, DataFactory_t>::Update(State& S, const Key& request) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Secondary evaluator \"" << my_key_ << "\" requested by "
               << request << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  for (auto& dep : dependencies_) {
    update |= S.GetEvaluator(dep.first, dep.second)->Update(S, Keys::getRequest(my_key_, my_tag_));
  }

  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Updating \"" << vo_.color("green") << my_key_ 
                 << vo_.reset() << "\" value... " << std::endl;
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
        *vo_.os() << "\"" << my_key_ << "\" has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "\"" << my_key_ << "\" has not changed... " << std::endl;
      }
      return false;
    }
  }
}


// ---------------------------------------------------------------------------
// Answers the question, has this data's derivative with respect to Key
// wrt_key changed since it was last requested by a requestor.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
bool EvaluatorSecondary<Data_t, DataFactory_t>::UpdateDerivative(State& S,
        const Key& requestor, const Key& wrt_key, const Key& wrt_tag) {
  ASSERT(IsDependency(S, wrt_key, wrt_tag));

  Teuchos::OSTab tab = vo_.getOSTab();
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Secondary Variable " << my_key_ << ":" << my_tag_ << " derivative with respect to "
              << wrt_key << ":" << wrt_tag << " requested by " << requestor;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;

  // -- must update if our our dependencies have changed, as these affect the partial derivatives
  Key my_request = Key{"d"}+Keys::getRequest(my_key_, my_tag_)+"_d"+Keys::getRequest(wrt_key, wrt_tag);
  update |= Update(S, my_request);

  // -- must update if any of our dependencies' derivatives have changed
  for (auto& dep : dependencies_) {
    if (S.GetEvaluator(dep.first, dep.second)->IsDependency(S, wrt_key, wrt_tag)) {
      update |= S.GetEvaluator(dep.first, dep.second)->UpdateDerivative(S, my_request, wrt_key, wrt_tag);
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
template<typename Data_t, typename DataFactory_t>
void EvaluatorSecondary<Data_t, DataFactory_t>::Update_(State& S) {
  // pull my variables out of state
  Data_t& myfield = S.GetW<Data_t>(my_key_, my_tag_, my_key_);

  // call the evaluate method
  Evaluate_(S, myfield);
}


// ---------------------------------------------------------------------------
// Updates the derivative
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
void EvaluatorSecondary<Data_t, DataFactory_t>::UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) {
  Key value = Keys::getRequest(my_key_, my_tag_);
  Key wrt = Keys::getRequest(wrt_key, wrt_tag);
  Errors::Message message;
  message << "Derivative requested for d" << value << "_d" << wrt
          << ", but derivatives are not implemented for this type.  Likely a missed template specialization?";
  throw(message);
}

// ---------------------------------------------------------------------------
// Does numerical FD to check the derivative
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
void EvaluatorSecondary<Data_t, DataFactory_t>::CheckDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) {
  Key value = Keys::getRequest(my_key_, my_tag_);
  Key wrt = Keys::getRequest(wrt_key, wrt_tag);
  Errors::Message message;
  message << "Derivative check requested for d" << value << "_d" << wrt
          << ", but FD is not implemented for this type.  Likely a missed template specialization or a lazy programmer?";
  throw(message);
}


// ---------------------------------------------------------------------------
// Does this evaluator depend upon the variable key@tag?
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
bool
EvaluatorSecondary<Data_t, DataFactory_t>::IsDependency(const State& S,
        const Key& key, const Key& tag) const {
  if (std::find(dependencies_.begin(), dependencies_.end(), std::make_pair(key,tag)) != dependencies_.end() ) {
    return true;
  } else {
    for (auto& dep : dependencies_) {
      if (S.GetEvaluator(dep.first, dep.second)->IsDependency(S,key,tag)) {
        return true;
      }
    }
  }
  return false;
}


// ---------------------------------------------------------------------------
// Does this evaluator provide the requested key@tag?
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
bool
EvaluatorSecondary<Data_t, DataFactory_t>::ProvidesKey(const Key& key, const Key& tag) const {
  return key == my_key_ && tag == my_tag_;
}


// ---------------------------------------------------------------------------
// Tracks through the dag ensuring ensuring requirements are met.
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
void
EvaluatorSecondary<Data_t, DataFactory_t>::EnsureCompatibility(State& S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(my_key_ != std::string(""));
  auto my_fac = S.Require<Data_t,DataFactory_t>(my_key_, my_tag_, my_key_);
  for (auto& dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);
  
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);
}



// ---------------------------------------------------------------------------
// Debugging
// ---------------------------------------------------------------------------
template<typename Data_t, typename DataFactory_t>
std::string
EvaluatorSecondary<Data_t, DataFactory_t>::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: secondary" << std::endl;
  for (auto& dep : dependencies_) {
    result << "  Dep: " << dep.first << "," << dep.second << std::endl;
  }
  result << std::endl;
  return result.str();
}




