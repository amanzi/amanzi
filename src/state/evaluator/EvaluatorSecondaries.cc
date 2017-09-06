/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default field evaluator base class.  A Evaluator is a node in the dependency
graph.

------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "EvaluatorSecondaries.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
EvaluatorSecondaries::EvaluatorSecondaries(
            Teuchos::ParameterList& plist) :
    my_tag_(plist.get<std::string>("tag","")),
    vo_(Keys::cleanPListName(plist.name()), plist),
    plist_(plist),
    check_derivative_(plist.get<bool>("check derivatives", false))
{
  // process the plist
  if (plist_.isParameter("evaluator names")) {
    Teuchos::Array<std::string> names =
        plist_.get<Teuchos::Array<std::string> >("evaluator names");
    for (Teuchos::Array<std::string>::iterator name=names.begin();
         name!=names.end(); ++name) {
      my_keys_.push_back(*name);
    }
  }

  if (plist_.isParameter("evaluator dependencies")) {
    Teuchos::Array<std::string> deps =
        plist_.get<Teuchos::Array<std::string> >("evaluator dependencies");
    for (Teuchos::Array<std::string>::iterator dep=deps.begin();
         dep!=deps.end(); ++dep) {
      dependencies_.insert(*dep);
    }
  }
}


Evaluator& EvaluatorSecondaries::operator=(const Evaluator& other) {
  if (this != &other) {
    const EvaluatorSecondaries* other_p =
        dynamic_cast<const EvaluatorSecondaries*>(&other);
    ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

EvaluatorSecondaries&
EvaluatorSecondaries::operator=(const EvaluatorSecondaries& other) {
  if (this != &other) {
    ASSERT(my_keys_[0] == other.my_keys_[0]);
    requests_ = other.requests_;
    deriv_requests_ = other.deriv_requests_;
  }
  return *this;
}


// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool EvaluatorSecondaries::Update(State& S,
        const Key& request) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Secondary Variable " << my_keys_[0] << " requested by "
          << request << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S.GetEvaluator(*dep)->Update(S, *my_keys_.begin());
  }

  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Updating " << my_keys_[0] << " value... " << std::endl;
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
        *vo_.os() << my_keys_[0]  << " has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << my_keys_[0] << " has not changed... " << std::endl;
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
bool EvaluatorSecondaries::UpdateDerivative(
    State& S, const Key& request, const Key& wrt_key) {
  Teuchos::OSTab tab = vo_.getOSTab();

  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    *vo_.os() << "Secondary Variable d" << my_keys_[0] << "_d" << wrt_key
          << " requested by " << request << std::endl;
  }

  // If wrt_key is not a dependency, no need to differentiate.
  if (!IsDependency(S, wrt_key)) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << wrt_key << " is not a dependency... " << std::endl;
    }
    return false;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  // -- must update if our our dependencies have changed, as these affect the partial derivatives
  update |= Update(S, my_keys_[0]);

  // -- must update if any of our dependencies derivatives have changed
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S.GetEvaluator(*dep)->UpdateDerivative(S, *my_keys_.begin(), wrt_key);
  }

  // Do the update
  std::pair<Key,Key> deriv_request(request, wrt_key);
  if (update) {
    if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
      *vo_.os() << "Updating d" << my_keys_[0] << "_d" << wrt_key << "... " << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateFieldDerivative_(S, wrt_key);
    deriv_requests_.clear();
    deriv_requests_.insert(deriv_request);
    return true;
  } else {
    // Otherwise, see if we have filled this request already.
    if (deriv_requests_.find(deriv_request) == deriv_requests_.end()) {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "d" << my_keys_[0] << "_d" << wrt_key << " has changed, but no need to update... " << std::endl;
      }

      deriv_requests_.insert(deriv_request);
      return true;
    } else {
      if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
        *vo_.os() << "d" << my_keys_[0] << "_d" << wrt_key << " has not changed... " << std::endl;
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
  for (std::vector<Key>::const_iterator key=my_keys_.begin(); key!=my_keys_.end(); ++key) {
    myfields.push_back(S.GetPtrW<CompositeVector>(*key, my_tag_, *key).ptr());
  }

  // call the evaluate method
  EvaluateField_(S, myfields);
}


// ---------------------------------------------------------------------------
// Updates the total derivative d(my_field)/d(wrt_field) in state S.
//
// (using p to indicate partial derivative...)
// d(my_field)/d(wrt_field) = p(my_field)/p(wrt_field) + sum_(dep) d(dep)/d(wrt_field)
// ---------------------------------------------------------------------------
void EvaluatorSecondaries::UpdateFieldDerivative_(State& S,
        const Key& wrt_key) {
  std::vector<Teuchos::Ptr<CompositeVector> > dmys;
  for (auto my_key=my_keys_.begin(); my_key!=my_keys_.end(); ++my_key) {
    Key dmy_key = std::string("d")+*my_key+std::string("_d")+wrt_key;
    Teuchos::RCP<CompositeVector> dmy;
    if (S.HasData(dmy_key, my_tag_)) {
      // Get the field...
      auto dmy = S.GetPtrW<CompositeVector>(dmy_key, my_tag_, *my_key);
      dmy->PutScalar(0.0);
      dmys.push_back(dmy.ptr());
    } else {
      // or create the field.  Note we have to do extra work that is normally
      // done by State in initialize.
      auto my_fac =
          S.Require<CompositeVector,CompositeVectorSpace>(*my_key, my_tag_);
      auto new_fac =
          S.Require<CompositeVector,CompositeVectorSpace>(dmy_key, my_tag_, *my_key);
      new_fac.Update(my_fac);
      dmy = Teuchos::rcp(new CompositeVector(new_fac));
      S.SetPtr(dmy_key, my_tag_, *my_key, dmy);
      S.GetRecordW(dmy_key, my_tag_, *my_key).set_initialized();
      S.GetRecordW(dmy_key, my_tag_, *my_key).set_io_vis(false);
      S.GetRecordW(dmy_key, my_tag_, *my_key).set_io_checkpoint(false);

      dmy->PutScalar(0.0);
      dmys.push_back(dmy.ptr());
    }
  }

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {

    // -- allocate a tmp set of vectors for dF_i/ddep
    std::vector<Teuchos::Ptr<CompositeVector> > tmp(dmys.size());
    for (int i=0; i!=dmys.size(); ++i) {
      tmp[i] = Teuchos::ptr(new CompositeVector(*dmys[i]));
    }

    if (wrt_key == *dep) {
      EvaluateFieldPartialDerivative_(S, wrt_key, tmp);
      for (int i=0; i!=dmys.size(); ++i) {
        // partial F_i / partial x
        dmys[i]->Update(1.0, *tmp[i], 1.0);
      }
    } else if (S.GetEvaluator(*dep)->IsDependency(S, wrt_key)) {
      EvaluateFieldPartialDerivative_(S, *dep, tmp);
      for (int i=0; i!=dmys.size(); ++i) {
        // partial F_i / partial dep * ddep/dx
        // -- ddep/dx
        Key ddep_key = std::string("d")+*dep+std::string("_d")+wrt_key;
        auto ddep = S.Get<CompositeVector>(ddep_key, my_tag_);
        dmys[i]->Multiply(1.0, ddep, *tmp[i], 1.0);
      }
    }

    // clean up tmp as it goes out of scope
    for (int i=0; i!=dmys.size(); ++i) {
      delete tmp[i].get();
    }
  }
  if (check_derivative_) {
    CheckDerivative_(S, wrt_key);
  }
}


inline
bool EvaluatorSecondaries::IsDependency(const State& S,
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
bool EvaluatorSecondaries::ProvidesKey(const Key& key) const {
  return (std::find(my_keys_.begin(), my_keys_.end(), key) != my_keys_.end());
}


void EvaluatorSecondaries::EnsureCompatibility(State& S) {
  CompositeVectorSpace* master_fac = nullptr;
  for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
       my_key!=my_keys_.end(); ++my_key) {
    // Ensure my field exists, and claim ownership.
    auto my_fac =
        S.Require<CompositeVector,CompositeVectorSpace>(*my_key, my_tag_, *my_key);
    S.GetRecordW(*my_key, my_tag_, *my_key).set_initialized();


    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>(std::string("visualize ")+*my_key, true);
    S.GetRecordW(*my_key, my_tag_, *my_key).set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+*my_key, false);
    S.GetRecordW(*my_key, my_tag_, *my_key).set_io_checkpoint(checkpoint_my_key);

    // Select a "master factory" to ensure commonality of all of my factories.
    if (my_fac.Mesh() != Teuchos::null) {
      if (master_fac == nullptr) {
        // no master yet, take this one
        master_fac = &my_fac;
      } else {
        // existing master, ensure it includes these requirements
        master_fac->Update(my_fac);
      }
    }
  }

  if (master_fac == nullptr) {
    // No requirements have been set, so we'll have to hope they get set by an
    // evaluator that depends upon this evaluator.
  } else {
    // Create an unowned factory to check my dependencies.
    CompositeVectorSpace dep_fac(*master_fac);
    dep_fac.SetOwned(false);

    // Loop over my keys, ensuring they meet the requirements of the master.
    for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
         my_key!=my_keys_.end(); ++my_key) {
      auto my_fac =
          S.Require<CompositeVector,CompositeVectorSpace>(*my_key, my_tag_);
      my_fac.Update(dep_fac);
    }

    // Loop over my dependencies, ensuring they meet my requirements.
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

void EvaluatorSecondaries::CheckDerivative_(State& S,
        const Key& wrt_key) {
  double eps = 1.e-10;

  // get unperturbed dependency and set up perturbed dependency
  auto dep_0 = S.GetPtrW<CompositeVector>(wrt_key, my_tag_,
          S.GetRecord(wrt_key, my_tag_).owner());
  auto dep_1 = Teuchos::rcp(new CompositeVector(*dep_0));
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
  std::vector<Teuchos::RCP<CompositeVector> > mys_0;
  std::vector<Teuchos::RCP<CompositeVector> > mys_1;

  for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
       my_key!=my_keys_.end(); ++my_key) {
    auto my_0 = S.GetPtrW<CompositeVector>(*my_key, my_tag_, *my_key);
    auto my_1 = Teuchos::rcp(new CompositeVector(*my_0));
    S.SetPtr(*my_key, my_tag_, *my_key, my_1);
    mys_0.push_back(my_0);
    mys_1.push_back(my_1);
  }
  Update(S, "test");

  for (int i=0; i!=my_keys_.size(); ++i) {
    Key my_key = my_keys_[i];
    Key dmy_key = std::string("d")+my_key+std::string("_d")+wrt_key;
    auto dmy = S.Get<CompositeVector>(dmy_key);

    // calculate differences in derivative estimates
    std::cout.precision(10);
    std::cout << "ERROR (eps=" << shift << ") for d_" << my_key << "_d" << wrt_key << std::endl;
    std::cout << "  dmy: " << dmy("cell",0) << ", " << dmy("cell",99) << std::endl;
    std::cout << "  my0 (" << (*dep_0)("cell",0) << ") = " << (*mys_0[i])("cell",0) << ", "
              << "my0 (" << (*dep_0)("cell",99) << ") = " << (*mys_0[i])("cell",99) << std::endl;
    std::cout << "  my1 (" << (*dep_1)("cell",0) << ") = " << (*mys_1[i])("cell",0) << ", "
              << "my1 (" << (*dep_1)("cell",99) << ") = " << (*mys_1[i])("cell",99) << std::endl;
    mys_1[i]->Update(-1.0/shift, *mys_0[i], 1.0/shift);
    std::cout << "  (my1-my0)/(dep1-dep0) = " << (*mys_1[i])("cell", 0) << ", = " << (*mys_1[i])("cell",99) << std::endl;

    // calculate relative error
    mys_1[i]->Update(-1.0, dmy, 1.0);
    std::cout << "  error (abs) = " << (*mys_1[i])("cell", 0) << ", = " << (*mys_1[i])("cell",99) << std::endl;

    double error1, error2;
    mys_1[i]->NormInf(&error1);
    dmy.NormInf(&error2);
    double error = error1 / (error2 + shift);
    std::cout << "ERROR (eps=" << shift << ") for d_" << my_key << "_d" << wrt_key << ":  " << error << std::endl;

    // revert pointers in state to old data
    S.SetPtr(my_key, my_tag_, my_key, mys_0[i]);
    S.SetPtr(wrt_key, my_tag_, S.GetRecord(wrt_key, my_tag_).owner(), dep_0);
  }

  // revert and clean up
  pfe->SetChanged();
  Update(S, "test");
}

std::string
EvaluatorSecondaries::WriteToString() const {
  std::stringstream result;
  for (std::vector<Key>::const_iterator my=my_keys_.begin();
       my!=my_keys_.end(); ++my) {
    result << *my << ",";
  }
  result << std::endl
         << "  Type: secondary" << std::endl;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep != dependencies_.end(); ++dep) {
    result << "  Dep: " << *dep << std::endl;
  }
  result << std::endl;
  return result.str();
}
} // namespace
