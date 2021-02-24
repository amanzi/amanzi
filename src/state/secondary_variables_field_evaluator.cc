/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default field evaluator base class.  A FieldEvaluator is a node in the Phalanx-like
dependency tree.

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
SecondaryVariablesFieldEvaluator::SecondaryVariablesFieldEvaluator(
            Teuchos::ParameterList& plist) :
    FieldEvaluator(plist)
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

  check_derivative_ = plist_.get<bool>("check derivatives", false);

  nonlocal_dependencies_ = plist_.get<bool>("includes non-rank-local dependencies", false);
}

SecondaryVariablesFieldEvaluator::SecondaryVariablesFieldEvaluator(
        const SecondaryVariablesFieldEvaluator& other) :
    FieldEvaluator(other),
    my_keys_(other.my_keys_),
    dependencies_(other.dependencies_),
    check_derivative_(other.check_derivative_),
    nonlocal_dependencies_(other.nonlocal_dependencies_)
{}


void SecondaryVariablesFieldEvaluator::operator=(const FieldEvaluator& other) {
#if ENABLE_DBC
  const SecondaryVariablesFieldEvaluator* other_p =
    dynamic_cast<const SecondaryVariablesFieldEvaluator*>(&other);
  AMANZI_ASSERT(other_p != NULL);
#else
  const SecondaryVariablesFieldEvaluator* other_p =
    static_cast<const SecondaryVariablesFieldEvaluator*>(&other);
#endif
  AMANZI_ASSERT(my_keys_[0] == other_p->my_keys_[0]);

  requests_ = other_p->requests_;
  deriv_requests_ = other_p->deriv_requests_;
}


// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool SecondaryVariablesFieldEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Secondary Variable " << my_keys_[0] << " requested by "
          << request << std::endl;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    bool updated = S->GetFieldEvaluator(*dep)->HasFieldChanged(S, *my_keys_.begin());
    if (updated && vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  must update " << my_keys_[0] << " due to changes in " << *dep << std::endl;
    update |= updated;
  }

  // check if nonlocal for changes in offprocess dependencies
  if (nonlocal_dependencies_) {
    int update_l = update;
    int update_g = 0;
    S->GetFieldData(my_keys_[0])->Comm()->MaxAll(&update_l, &update_g, 1);
    update |= update_g;
  }
  
  if (update) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Updating " << my_keys_[0] << " value... " << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateField_(S);
    requests_.clear();
    requests_.insert(request);
    return true;
  } else {
    // Otherwise, see if we have filled this request already.
    if (requests_.find(request) == requests_.end()) {
      requests_.insert(request);
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << my_keys_[0]  << " has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << my_keys_[0] << " has not changed... " << std::endl;
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
bool SecondaryVariablesFieldEvaluator::HasFieldDerivativeChanged(
    const Teuchos::Ptr<State>& S, Key request, Key wrt_key) {
  AMANZI_ASSERT(!request.empty());
  AMANZI_ASSERT(!wrt_key.empty());

  Teuchos::OSTab tab = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Secondary Variable " << Keys::getDerivKey(my_keys_[0], wrt_key)
          << " requested by " << request << std::endl;
  }

  // If wrt_key is not a dependency, no need to differentiate.
  if (!IsDependency(S, wrt_key)) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << wrt_key << " is not a dependency... " << std::endl;
    }
    return false;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  // -- must update if our our dependencies have changed, as these affect the partial derivatives
  update |= HasFieldChanged(S, my_keys_[0]);

  // -- must update if any of our dependencies derivatives have changed
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S->GetFieldEvaluator(*dep)->HasFieldDerivativeChanged(S, *my_keys_.begin(), wrt_key);
  }
  // check if nonlocal for changes in offprocess dependencies
  if (nonlocal_dependencies_) {
    int update_l = update;
    int update_g = 0;
    S->GetFieldData(my_keys_[0])->Comm()->MaxAll(&update_l, &update_g, 1);
    update |= update_g;
  }

  // Do the update
  std::pair<Key,Key> deriv_request(request, wrt_key);
  if (update) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Updating " << Keys::getDerivKey(my_keys_[0], wrt_key) << "... " << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateFieldDerivative_(S, wrt_key);
    deriv_requests_.clear();
    deriv_requests_.insert(deriv_request);
    return true;
  } else {
    // Otherwise, see if we have filled this request already.
    if (deriv_requests_.find(deriv_request) == deriv_requests_.end()) {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << Keys::getDerivKey(my_keys_[0], wrt_key) << " has changed, but no need to update... " << std::endl;
      }

      deriv_requests_.insert(deriv_request);
      return true;
    } else {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << Keys::getDerivKey(my_keys_[0], wrt_key) << " has not changed... " << std::endl;
      }
      return false;
    }
  }
}


// ---------------------------------------------------------------------------
// Updates the field value in state S.
// ---------------------------------------------------------------------------
void SecondaryVariablesFieldEvaluator::UpdateField_(const Teuchos::Ptr<State>& S) {
  // pull my variables out of state
  std::vector<Teuchos::Ptr<CompositeVector> > myfields;
  for (std::vector<Key>::const_iterator key=my_keys_.begin(); key!=my_keys_.end(); ++key) {
    myfields.push_back(S->GetFieldData(*key, *key).ptr());
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
void SecondaryVariablesFieldEvaluator::UpdateFieldDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key) {
  std::vector<Teuchos::Ptr<CompositeVector> > dmys;
  for (std::vector<Key>::const_iterator my_key=my_keys_.begin(); my_key!=my_keys_.end(); ++my_key) {
    Key dmy_key = Keys::getDerivKey(*my_key, wrt_key);
    if (S->HasField(dmy_key)) {
      // Get the field...
      Teuchos::RCP<CompositeVector> dmy = S->GetFieldData(dmy_key, *my_key);
      dmy->PutScalar(0.0);
      dmys.push_back(dmy.ptr());
    } else {
      // or create the field.  Note we have to do extra work that is normally
      // done by State in initialize.
      Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(*my_key);
      Teuchos::RCP<CompositeVectorSpace> new_fac =
        S->RequireField(dmy_key, *my_key);
      new_fac->Update(*my_fac);
      auto dmy = Teuchos::rcp(new CompositeVector(*new_fac));
      S->SetData(dmy_key, *my_key, dmy);
      S->GetField(dmy_key,*my_key)->set_initialized();
      S->GetField(dmy_key,*my_key)->set_io_vis(plist_.get<bool>("visualize derivative", false));
      S->GetField(dmy_key,*my_key)->set_io_checkpoint(plist_.get<bool>("checkpoint derivative", false));

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
    } else if (S->GetFieldEvaluator(*dep)->IsDependency(S, wrt_key)) {
      EvaluateFieldPartialDerivative_(S, *dep, tmp);
      for (int i=0; i!=dmys.size(); ++i) {
        // partial F_i / partial dep * ddep/dx
        // -- ddep/dx
        Key ddep_key = Keys::getDerivKey(*dep, wrt_key);
        Teuchos::RCP<const CompositeVector> ddep = S->GetFieldData(ddep_key);

        dmys[i]->Multiply(1.0, *ddep, *tmp[i], 1.0);
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


bool SecondaryVariablesFieldEvaluator::IsDependency(const Teuchos::Ptr<State>& S,
        Key key) const {
  if (dependencies_.find(key) != dependencies_.end()) {
    return true;
  } else {
    for (KeySet::const_iterator dep=dependencies_.begin();
         dep!=dependencies_.end(); ++dep) {
      if (S->GetFieldEvaluator(*dep)->IsDependency(S,key)) {
        return true;
      }
    }
  }
  return false;
}


bool SecondaryVariablesFieldEvaluator::ProvidesKey(Key key) const {
  return (std::find(my_keys_.begin(), my_keys_.end(), key) != my_keys_.end());
}


void SecondaryVariablesFieldEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<CompositeVectorSpace> master_fac;
  for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
       my_key!=my_keys_.end(); ++my_key) {
    // Ensure my field exists, and claim ownership.
    Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(*my_key, *my_key);

    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>("visualize", true);
    S->GetField(*my_key, *my_key)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
    S->GetField(*my_key, *my_key)->set_io_checkpoint(checkpoint_my_key);

    // Select a "master factory" to ensure commonality of all of my factories.
    if (my_fac->Mesh() != Teuchos::null) {
      if (master_fac == Teuchos::null) {
        // no master yet, take this one
        master_fac = my_fac;
      } else {
        // existing master, ensure it includes these requirements
        master_fac->Update(*my_fac);
      }
    }
  }

  if (master_fac == Teuchos::null) {
    // No requirements have been set, so we'll have to hope they get set by an
    // evaluator that depends upon this evaluator.
  } else {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace(*master_fac));
    dep_fac->SetOwned(false);

    // Loop over my keys, ensuring they meet the requirements of the master.
    for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
         my_key!=my_keys_.end(); ++my_key) {
      Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(*my_key);
      my_fac->Update(*dep_fac);
    }

    // Loop over my dependencies, ensuring they meet my requirements.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(*key);
      fac->Update(*dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
    }
  }
}

void SecondaryVariablesFieldEvaluator::CheckDerivative_(const Teuchos::Ptr<State>& S,
                                                       Key wrt_key) {
  double eps = 1.e-10;

  // get unperturbed dependency and set up perturbed dependency
  Teuchos::RCP<CompositeVector> dep_0 = S->GetFieldData(wrt_key, S->GetField(wrt_key)->owner());
  Teuchos::RCP<CompositeVector> dep_1 = Teuchos::rcp(new CompositeVector(*dep_0));
  *dep_1 = *dep_0;
  double norm;
  dep_0->Norm2(&norm);
  double shift = (eps*abs(norm) + eps);
  dep_1->Shift(shift);

  S->SetData(wrt_key, S->GetField(wrt_key)->owner(), dep_1);
  Teuchos::RCP<FieldEvaluator> fe = S->GetFieldEvaluator(wrt_key);
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
  pfe->SetFieldAsChanged(S);

  // get unperturbed value and calulate perturbed value
  std::vector<Teuchos::RCP<CompositeVector> > mys_0;
  std::vector<Teuchos::RCP<CompositeVector> > mys_1;

  for (std::vector<Key>::const_iterator my_key=my_keys_.begin();
       my_key!=my_keys_.end(); ++my_key) {
    Teuchos::RCP<CompositeVector> my_0 = S->GetFieldData(*my_key, *my_key);
    Teuchos::RCP<CompositeVector> my_1 = Teuchos::rcp(new CompositeVector(*my_0));
    S->SetData(*my_key, *my_key, my_1);
    mys_0.push_back(my_0);
    mys_1.push_back(my_1);
  }
  HasFieldChanged(S, "test");

  for (int i=0; i!=my_keys_.size(); ++i) {
    Key my_key = my_keys_[i];
    Key dmy_key = Keys::getDerivKey(my_key, wrt_key);
    Teuchos::RCP<const CompositeVector> dmy = S->GetFieldData(dmy_key);

    // calculate differences in derivative estimates
    std::cout.precision(10);
    std::cout << "ERROR (eps=" << shift << ") for " << Keys::getDerivKey(my_key, wrt_key) << std::endl;
    std::cout << "  dmy: " << (*dmy)("cell",0) << ", " << (*dmy)("cell",99) << std::endl;
    std::cout << "  my0 (" << (*dep_0)("cell",0) << ") = " << (*mys_0[i])("cell",0) << ", "
              << "my0 (" << (*dep_0)("cell",99) << ") = " << (*mys_0[i])("cell",99) << std::endl;
    std::cout << "  my1 (" << (*dep_1)("cell",0) << ") = " << (*mys_1[i])("cell",0) << ", "
              << "my1 (" << (*dep_1)("cell",99) << ") = " << (*mys_1[i])("cell",99) << std::endl;
    mys_1[i]->Update(-1.0/shift, *mys_0[i], 1.0/shift);
    std::cout << "  (my1-my0)/(dep1-dep0) = " << (*mys_1[i])("cell", 0) << ", = " << (*mys_1[i])("cell",99) << std::endl;

    // calculate relative error
    mys_1[i]->Update(-1.0, *dmy, 1.0);
    std::cout << "  error (abs) = " << (*mys_1[i])("cell", 0) << ", = " << (*mys_1[i])("cell",99) << std::endl;

    double error1, error2;
    mys_1[i]->NormInf(&error1);
    dmy->NormInf(&error2);
    double error = error1 / (error2 + shift);
    std::cout << "ERROR (eps=" << shift << ") for " << Keys::getDerivKey(my_key, wrt_key) << ":  " << error << std::endl;

    // revert pointers in state to old data
    S->SetData(my_key, my_key, mys_0[i]);
    S->SetData(wrt_key, S->GetField(wrt_key)->owner(), dep_0);
  }

  // revert and clean up
  pfe->SetFieldAsChanged(S);
  HasFieldChanged(S, "test");
}

std::string
SecondaryVariablesFieldEvaluator::WriteToString() const {
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
