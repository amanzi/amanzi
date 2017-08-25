/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default field evaluator base class.  A Evaluator is a node in the dependency
graph.

------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
EvaluatorSecondary::EvaluatorSecondary(Teuchos::ParameterList& plist) :
    my_key_(plist.name()),
    my_tag_(plist.get<std::string>("tag","")),
    vo_(plist.name(), plist),
    plist_(plist),
    check_derivative_(plist.get<bool>("check derivatives", false))
{
  if (plist_.isParameter("evaluator dependencies")) {
    Teuchos::Array<std::string> deps =
        plist_.get<Teuchos::Array<std::string> >("evaluator dependencies");
    for (auto dep : deps) {
      dependencies_.insert(dep);
    }
  }
}


// -----------------------------------------------------------------------------
// Assignment operators
// -----------------------------------------------------------------------------
Evaluator& EvaluatorSecondary::operator=(const Evaluator& other) {
  if (this != &other) {
    const EvaluatorSecondary* other_p =
        dynamic_cast<const EvaluatorSecondary*>(&other);
    ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


// -----------------------------------------------------------------------------
// Assignment operators
// -----------------------------------------------------------------------------
EvaluatorSecondary& EvaluatorSecondary::operator=(const EvaluatorSecondary& other) {
  if (this != &other) {
    ASSERT(my_key_ == other.my_key_);
    requests_ = other_p->requests_;
    deriv_requests_ = other_p->deriv_requests_;
  }
  return *this;
}

// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool EvaluatorSecondary::Update(State& S,
        const Key& request) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Secondary field \"" << my_key_ << "\" requested by "
               << request << std::endl;
  }


  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S.GetEvaluator(*dep, my_tag_)->Update(S, my_key_);
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
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
bool EvaluatorSecondary::UpdateDerivative(State& S,
        const Key& request, const Key& wrt_key) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Secondary Variable d" << my_key_ << "_d" << wrt_key
          << " requested by " << request << std::endl;
  }

  // If wrt_key is not a dependency, no need to differentiate.
  if (!IsDependency(S, wrt_key)) {
    return false;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  // -- must update if our our dependencies have changed, as these affect the partial derivatives
  update |= Update(S, my_key_);

  // -- must update if any of our dependencies derivatives have changed
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S.GetEvaluator(*dep, my_tag_)->UpdateDerivative(S, my_key_, wrt_key);
  }

  // Do the update
  std::pair<Key,Key> deriv_request(request, wrt_key);
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Updating d" << my_key_ << "_d" << wrt_key << "... " << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateFieldDerivative_(S, wrt_key);
    deriv_requests_.clear();
    deriv_requests_.insert(deriv_request);
    return true;
  } else {
    // Otherwise, simply service the request
    if (deriv_requests_.find(deriv_request) == deriv_requests_.end()) {
      deriv_requests_.insert(deriv_request);
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "d" << my_key_ << "_d" << wrt_key << " has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "d" << my_key_ << "_d" << wrt_key << " has not changed... " << std::endl;
      }
      return false;
    }
  }
}


// ---------------------------------------------------------------------------
// Updates the field value in state S.
// ---------------------------------------------------------------------------
void EvaluatorSecondary::Update_(State& S) {
  // pull my variables out of state
  Teuchos::RCP<CompositeVector> myfield =
      S.GetPtrW<CompositeVector>(my_key_, my_tag_, my_key_);

  // call the evaluate method
  EvaluateField_(S, myfield.ptr());
}


// ---------------------------------------------------------------------------
// Updates the total derivative d(my_field)/d(wrt_field) in state S.
//
// (using p to indicate partial derivative...)
// d(my_field)/d(wrt_field) = p(my_field)/p(wrt_field) + sum_(dep) d(dep)/d(wrt_field)
// ---------------------------------------------------------------------------
void EvaluatorSecondary::UpdateFieldDerivative_(State& S,
        const Key& wrt_key) {
  Key dmy_key = std::string("d")+my_key_+std::string("_d")+wrt_key;
  Teuchos::RCP<CompositeVector> dmy;
  if (S.HasData(dmy_key, my_tag_)) {
    // Get the field...
    dmy = S.GetPtrW<CompositeVector>(dmy_key, my_tag_, my_key_);
  } else {
    // or create the field.  Note we have to do extra work that is normally
    // done by State in initialize.
    auto my_fac =
        S.Require<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_);
    auto new_fac =
        S.Require<CompositeVector,CompositeVectorSpace>(dmy_key, my_tag_, my_key_);
    new_fac.Update(my_fac);
    new_fac.SetGhosted(my_fac.Ghosted());
    dmy = Teuchos::rcp(new CompositeVector(new_fac));
    S.SetPtr(dmy_key, my_tag_, my_key_, dmy);
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_initialized();
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_io_vis(false);
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_io_checkpoint(false);
  }

  // Now update the values.
  dmy->PutScalar(0.0);

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    Teuchos::RCP<CompositeVector> tmp = Teuchos::rcp(new CompositeVector(*dmy));
    if (wrt_key == *dep) {
      // partial F / partial x
      EvaluateFieldPartialDerivative_(S, wrt_key, tmp.ptr());
      dmy->Update(1.0, *tmp, 1.0);

      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << dmy_key << " = " << (*tmp)("cell",0) << std::endl;
      }

    } else if (S.GetEvaluator(*dep)->IsDependency(S, wrt_key)) {
      // partial F / partial dep * ddep/dx
      // -- ddep/dx
      Key ddep_key = std::string("d")+*dep+std::string("_d")+wrt_key;
      Teuchos::RCP<const CompositeVector> ddep = S.GetFieldData(ddep_key);
      // -- partial F / partial dep
      EvaluateFieldPartialDerivative_(S, *dep, tmp.ptr());

      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << ddep_key << " = " << (*ddep)("cell",0) << ", ";
        *vo_.os() << "d"<< my_key_ << "_d" << *dep << " = " << (*tmp)("cell",0) << std::endl;
      }

      dmy->Multiply(1.0, *ddep, *tmp, 1.0);
    }

    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << dmy_key << " = " << (*dmy)("cell",0) << std::endl;
    }
  }

  if (check_derivative_) {
    CheckDerivative_(S, wrt_key);
  }
}


inline
bool EvaluatorSecondary::IsDependency(const State& S,
        const Key& key) const {
  if (dependencies_.find(key) != dependencies_.end()) {
    return true;
  } else {
    for (KeySet::const_iterator dep=dependencies_.begin();
         dep!=dependencies_.end(); ++dep) {
      if (S.GetEvaluator(*dep)->IsDependency(S,key)) {
        return true;
      }
    }
  }
  return false;
}


inline
bool EvaluatorSecondary::ProvidesKey(const Key& key) const {
  return (key == my_key_);
}


void EvaluatorSecondary::EnsureCompatibility(State& S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(my_key_ != std::string(""));
  auto my_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    CompositeVectorSpace dep_fac(my_fac);
    dep_fac.SetOwned(false);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      auto fac = S.Require<CompositeVector,CompositeVectorSpace>(*key, my_tag_);
      fac.Update(dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S.RequireEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}


void EvaluatorSecondary::CheckDerivative_(State& S,
        const Key& wrt_key) {
  double eps = 1.e-10;

  Key dmy_key = std::string("d")+my_key_+std::string("_d")+wrt_key;
  const CompositeVector& dmy = S.Get<CompositeVector>(dmy_key, my_tag_);

  // get unperturbed dependency and set up perturbed dependency
  auto dep_0 = S.GetPtrW<CompositeVector>(wrt_key, my_tag_,
          S.GetRecord(wrt_key, my_tag_).owner());
  Teuchos::RCP<CompositeVector> dep_1 = Teuchos::rcp(new CompositeVector(*dep_0));
  *dep_1 = *dep_0;
  double norm;
  dep_0->Norm2(&norm);
  double shift = (eps*abs(norm) + eps);
  dep_1->Shift(shift);
  S.SetPtr(wrt_key, my_tag_, S.GetRecord(wrt_key, my_tag_).owner(), dep_1);
  Teuchos::RCP<Evaluator> fe = S.GetEvaluator(wrt_key, my_tag_);
  Teuchos::RCP<EvaluatorPrimary> pfe = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(fe);
  pfe->SetChanged();

  // get unperturbed value and calulate perturbed value
  Teuchos::RCP<CompositeVector> my_0 = S.GetPtrW<CompositeVector>(my_key_, my_tag_, my_key_);
  Teuchos::RCP<CompositeVector> my_1 = Teuchos::rcp(new CompositeVector(*my_0));
  S.SetPtr(my_key_, my_tag_, my_key_, my_1);
  Update(S, "test");

  // calculate differences in derivative estimates
  std::cout.precision(10);
  std::cout << "ERROR (eps=" << shift << ") for " << dmy_key << std::endl;
  std::cout << "  dmy: " << dmy("cell",0) << ", " << dmy("cell",0) << std::endl;
  std::cout << "  my0 (" << (*dep_0)("cell",0) << ") = " << (*my_0)("cell",0) << ", "
            << "my0 (" << (*dep_0)("cell",0) << ") = " << (*my_0)("cell",0) << std::endl;
  std::cout << "  my1 (" << (*dep_1)("cell",0) << ") = " << (*my_1)("cell",0) << ", "
            << "my1 (" << (*dep_1)("cell",0) << ") = " << (*my_1)("cell",0) << std::endl;
  my_1->Update(-1.0/shift, *my_0, 1.0/shift);
  std::cout << "  (my1-my0)/(dep1-dep0) = " << (*my_1)("cell", 0) << ", = " << (*my_1)("cell",0) << std::endl;

  // calculate relative error
  my_1->Update(-1.0, dmy, 1.0);
  std::cout << "  error (abs) = " << (*my_1)("cell", 0) << ", = " << (*my_1)("cell",0) << std::endl;

  double error1, error2;
  my_1->NormInf(&error1);
  dmy.NormInf(&error2);
  double error = error1 / (error2 + shift);
  std::cout << "ERROR (eps=" << shift << ") for d" << dmy_key << ":  " << error << std::endl;

  // revert pointers in state to old data
  S.SetPtr(my_key_, my_tag_, my_key_, my_0);
  S.SetPtr(wrt_key, my_tag_, S.GetRecord(wrt_key, my_tag_).owner(), dep_0);
  pfe->SetChanged();
  Update(S, "test");
}


std::string
EvaluatorSecondary::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: secondary" << std::endl;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep != dependencies_.end(); ++dep) {
    result << "  Dep: " << *dep << std::endl;
  }
  result << std::endl;
  return result.str();
}

} // namespace
