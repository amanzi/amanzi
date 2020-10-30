/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default field evaluator base class.  A FieldEvaluator is a node in the Phalanx-like
dependency tree.

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------

SecondaryVariableFieldEvaluator::SecondaryVariableFieldEvaluator(
        Teuchos::ParameterList& plist) :
    FieldEvaluator(plist)
{
  // process the plist
  if (plist_.isParameter("evaluator name")) {
    my_key_ = plist_.get<std::string>("evaluator name");
  }

  if (plist_.isParameter("evaluator dependencies")) {
    Teuchos::Array<std::string> deps = plist_.get<Teuchos::Array<std::string> >("evaluator dependencies");
    for (const auto& dep : deps) {
      dependencies_.insert(dep);
    }
  }

  check_derivative_ = plist_.get<bool>("check derivatives", false);

  nonlocal_dependencies_ = plist_.get<bool>("includes non-rank-local dependencies", false);
}

SecondaryVariableFieldEvaluator::SecondaryVariableFieldEvaluator(
        const SecondaryVariableFieldEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_),
    dependencies_(other.dependencies_),
    check_derivative_(other.check_derivative_),
    nonlocal_dependencies_(other.nonlocal_dependencies_)
{}


void SecondaryVariableFieldEvaluator::operator=(const FieldEvaluator& other) {
#if ENABLE_DBC
  const SecondaryVariableFieldEvaluator* other_p =
    dynamic_cast<const SecondaryVariableFieldEvaluator*>(&other);
  AMANZI_ASSERT(other_p != NULL);
#else
  const SecondaryVariableFieldEvaluator* other_p =
    static_cast<const SecondaryVariableFieldEvaluator*>(&other);
#endif
  AMANZI_ASSERT(my_key_ == other_p->my_key_);

  requests_ = other_p->requests_;
  deriv_requests_ = other_p->deriv_requests_;
}

// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool SecondaryVariableFieldEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Secondary field \"" << my_key_ << "\" requested by "
               << request << std::endl;
  }


  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S->GetFieldEvaluator(*dep)->HasFieldChanged(S, my_key_);
  }
  // check if nonlocal for changes in offprocess dependencies
  if (nonlocal_dependencies_) {
    int update_l = update;
    int update_g = 0;
    S->GetFieldData(my_key_)->Comm()->MaxAll(&update_l, &update_g, 1);
    update |= update_g;
  }

  if (update) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Updating \"" << vo_->color("green") << my_key_ 
                 << vo_->reset() << "\" value... " << std::endl;
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
        *vo_->os() << "\"" << my_key_ << "\" has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "\"" << my_key_ << "\" has not changed... " << std::endl;
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
bool SecondaryVariableFieldEvaluator::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {
  Teuchos::OSTab tab = vo_->getOSTab();

  AMANZI_ASSERT(!wrt_key.empty());
  AMANZI_ASSERT(!request.empty());
  
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Secondary Variable d" << my_key_ << "_d" << wrt_key
          << " requested by " << request << std::endl;
  }

  // If wrt_key is not a dependency, no need to differentiate.
  if (!IsDependency(S, wrt_key)) {
    return false;
  }

  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;
  // -- must update if our our dependencies have changed, as these affect the partial derivatives
  update |= HasFieldChanged(S, my_key_);

  // -- must update if any of our dependencies derivatives have changed
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    update |= S->GetFieldEvaluator(*dep)->HasFieldDerivativeChanged(S, my_key_, wrt_key);
  }

  // check if nonlocal for changes in offprocess dependencies
  if (nonlocal_dependencies_) {
    int update_l = update;
    int update_g = 0;
    S->GetFieldData(my_key_)->Comm()->MaxAll(&update_l, &update_g, 1);
    update |= update_g;
  }
  
  // Do the update
  std::pair<Key,Key> deriv_request(request, wrt_key);
  if (update) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Updating d" << my_key_ << "_d" << wrt_key << "... " << std::endl;
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
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "d" << my_key_ << "_d" << wrt_key << " has changed, but no need to update... " << std::endl;
      }
      return true;
    } else {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << "d" << my_key_ << "_d" << wrt_key << " has not changed... " << std::endl;
      }
      return false;
    }
  }
}


// ---------------------------------------------------------------------------
// Updates the field value in state S.
// ---------------------------------------------------------------------------
void SecondaryVariableFieldEvaluator::UpdateField_(const Teuchos::Ptr<State>& S) {
  // pull my variables out of state
  Teuchos::RCP<CompositeVector> myfield = S->GetFieldData(my_key_, my_key_);

  // call the evaluate method
  EvaluateField_(S, myfield.ptr());
}


// ---------------------------------------------------------------------------
// Updates the total derivative d(my_field)/d(wrt_field) in state S.
//
// (using p to indicate partial derivative...)
// d(my_field)/d(wrt_field) = p(my_field)/p(wrt_field) + sum_(dep) d(dep)/d(wrt_field)
// ---------------------------------------------------------------------------
void SecondaryVariableFieldEvaluator::UpdateFieldDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key) {
  Key dmy_key = std::string("d")+my_key_+std::string("_d")+wrt_key;
  Teuchos::RCP<CompositeVector> dmy;
  if (S->HasField(dmy_key)) {
    // Get the field...
    dmy = S->GetFieldData(dmy_key, my_key_);
  } else {
    // or create the field.  Note we have to do extra work that is normally
    // done by State in initialize.
    Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_);
    Teuchos::RCP<CompositeVectorSpace> new_fac =
      S->RequireField(dmy_key, my_key_);
    new_fac->Update(*my_fac);
    new_fac->SetGhosted(my_fac->Ghosted());
    dmy = Teuchos::rcp(new CompositeVector(*new_fac));
    S->SetData(dmy_key, my_key_, dmy);
    S->GetField(dmy_key,my_key_)->set_initialized();
    S->GetField(dmy_key,my_key_)->set_io_vis(plist_.get<bool>("visualize derivative", false));
    S->GetField(dmy_key,my_key_)->set_io_checkpoint(plist_.get<bool>("checkpoint derivative", false));
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

      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << dmy_key << " = " << (*tmp)("cell",0) << std::endl;
      }

    } else if (S->GetFieldEvaluator(*dep)->IsDependency(S, wrt_key)) {
      // partial F / partial dep * ddep/dx
      // -- ddep/dx
      Key ddep_key = std::string("d")+*dep+std::string("_d")+wrt_key;
      Teuchos::RCP<const CompositeVector> ddep = S->GetFieldData(ddep_key);
      // -- partial F / partial dep
      EvaluateFieldPartialDerivative_(S, *dep, tmp.ptr());

      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        *vo_->os() << ddep_key << " = " << (*ddep)("cell",0) << ", ";
        *vo_->os() << "d"<< my_key_ << "_d" << *dep << " = " << (*tmp)("cell",0) << std::endl;
      }

      dmy->Multiply(1.0, *ddep, *tmp, 1.0);
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << dmy_key << " = " << (*dmy)("cell",0) << std::endl;
    }
  }

  if (check_derivative_) {
    CheckDerivative_(S, wrt_key);
  }
}


bool SecondaryVariableFieldEvaluator::IsDependency(const Teuchos::Ptr<State>& S,
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


bool SecondaryVariableFieldEvaluator::ProvidesKey(Key key) const {
  return (key == my_key_);
}


void SecondaryVariableFieldEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_key_ != std::string(""));
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace(*my_fac));
    dep_fac->SetOwned(false);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& key : dependencies_) {
      if (key == my_key_) {
        Errors::Message msg;
        msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
      fac->Update(*dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (const auto& key : dependencies_) {
      S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
    }
  }
}


void SecondaryVariableFieldEvaluator::CheckDerivative_(const Teuchos::Ptr<State>& S,
                                                       Key wrt_key) {
  double eps = 1.e-10;

  Key dmy_key = std::string("d")+my_key_+std::string("_d")+wrt_key;
  Teuchos::RCP<const CompositeVector> dmy = S->GetFieldData(dmy_key);

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
  Teuchos::RCP<CompositeVector> my_0 = S->GetFieldData(my_key_, my_key_);
  Teuchos::RCP<CompositeVector> my_1 = Teuchos::rcp(new CompositeVector(*my_0));
  S->SetData(my_key_, my_key_, my_1);
  HasFieldChanged(S, "test");

  // calculate differences in derivative estimates
  std::cout.precision(10);
  std::cout << "ERROR (eps=" << shift << ") for " << dmy_key << std::endl;
  std::cout << "  dmy: " << (*dmy)("cell",0) << ", " << (*dmy)("cell",0) << std::endl;
  std::cout << "  my0 (" << (*dep_0)("cell",0) << ") = " << (*my_0)("cell",0) << ", "
            << "my0 (" << (*dep_0)("cell",0) << ") = " << (*my_0)("cell",0) << std::endl;
  std::cout << "  my1 (" << (*dep_1)("cell",0) << ") = " << (*my_1)("cell",0) << ", "
            << "my1 (" << (*dep_1)("cell",0) << ") = " << (*my_1)("cell",0) << std::endl;
  my_1->Update(-1.0/shift, *my_0, 1.0/shift);
  std::cout << "  (my1-my0)/(dep1-dep0) = " << (*my_1)("cell", 0) << ", = " << (*my_1)("cell",0) << std::endl;

  // calculate relative error
  my_1->Update(-1.0, *dmy, 1.0);
  std::cout << "  error (abs) = " << (*my_1)("cell", 0) << ", = " << (*my_1)("cell",0) << std::endl;

  double error1, error2;
  my_1->NormInf(&error1);
  dmy->NormInf(&error2);
  double error = error1 / (error2 + shift);
  std::cout << "ERROR (eps=" << shift << ") for d" << dmy_key << ":  " << error << std::endl;

  // revert pointers in state to old data
  S->SetData(my_key_, my_key_, my_0);
  S->SetData(wrt_key, S->GetField(wrt_key)->owner(), dep_0);
  pfe->SetFieldAsChanged(S);
  HasFieldChanged(S, "test");
}


std::string
SecondaryVariableFieldEvaluator::WriteToString() const {
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
