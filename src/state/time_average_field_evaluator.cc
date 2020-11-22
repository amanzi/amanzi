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
#include "time_average_field_evaluator.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------

TimeAverageFieldEvaluator::TimeAverageFieldEvaluator(
        Teuchos::ParameterList& plist) :
    FieldEvaluator(plist)
{
  // process the plist
  if (plist_.isParameter("evaluator name")) {
    my_key_ = plist_.get<std::string>("evaluator name");
  }


  base_key_ = plist_.get<std::string>("base field");
  dependencies_.insert(base_key_);


  current_bin_ = 0;
  number_bins_ = plist_.get<int>("number of bins");
  AMANZI_ASSERT(number_bins_ > 0);
  bin_period_ = plist_.get<double>("bin period");
  AMANZI_ASSERT(bin_period_ > 0);
  bin_time_.resize(number_bins_);
  update_time_ = 0.;

}

TimeAverageFieldEvaluator::TimeAverageFieldEvaluator(
        const TimeAverageFieldEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_),
    dependencies_(other.dependencies_),
    requests_(other.requests_),
    current_bin_(other.current_bin_),
    number_bins_(other.number_bins_)
{
  for (int i=0;i<bin_time_.size(); ++i) bin_time_[i] = other.bin_time_[i];
  field_time_bins_ = other.field_time_bins_;
  bin_period_ = other.bin_period_;
  update_time_ = other.update_time_;
  base_key_ = other.base_key_;
  dependencies_.insert(base_key_);
}


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------  
Teuchos::RCP<FieldEvaluator> TimeAverageFieldEvaluator::Clone() const {
  return Teuchos::rcp(new TimeAverageFieldEvaluator(*this));
}

void TimeAverageFieldEvaluator::operator=(const FieldEvaluator& other) {
#if ENABLE_DBC
  const TimeAverageFieldEvaluator* other_p =
    dynamic_cast<const TimeAverageFieldEvaluator*>(&other);
  AMANZI_ASSERT(other_p != NULL);
#else
  const TimeAverageFieldEvaluator* other_p =
    static_cast<const TimeAverageFieldEvaluator*>(&other);
#endif
  AMANZI_ASSERT(my_key_ == other_p->my_key_);

  current_bin_ = other_p->current_bin_;
  number_bins_ = other_p->number_bins_;  
  requests_ = other_p->requests_;

  bin_time_.resize(other_p->bin_time_.size());
  AMANZI_ASSERT(bin_time_.size() == other_p->bin_time_.size());
  if (bin_time_.size() > 0) {
    for (int i=0;i<bin_time_.size(); ++i) bin_time_[i] = other_p->bin_time_[i];
  }
  field_time_bins_ = other_p->field_time_bins_;
  bin_period_ = other_p->bin_period_;
  update_time_ = other_p->update_time_;  
}

// -----------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// -----------------------------------------------------------------------------
bool TimeAverageFieldEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S,
        Key request) {
  Teuchos::OSTab tab = vo_->getOSTab();

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "Time averaged field \"" << my_key_ << "\" requested by "
               << request << std::endl;
  }


  // Check if we need to update ourselves, and potentially update our dependencies.
  bool update = false;

  if (S->time() - update_time_ > 1e-10) update = true;
  
  if (update) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Updating \"" << vo_->color("green") << my_key_ 
                 << vo_->reset() << "\" value... " << std::endl;
    }

    // If so, update ourselves, empty our list of filled requests, and return.
    UpdateField_(S);
    requests_.clear();
    requests_.insert(request);
    update_time_ = S->time();
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
// Updates the field value in state S.
// ---------------------------------------------------------------------------
void TimeAverageFieldEvaluator::UpdateField_(const Teuchos::Ptr<State>& S) {
  // pull my variables out of state
  Teuchos::RCP<CompositeVector> myfield = S->GetFieldData(my_key_, my_key_);

  // call the evaluate method
  EvaluateField_(S, myfield.ptr());
}


void TimeAverageFieldEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                      const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> base_field= S->GetFieldData(base_key_);
  double norm[20], res_norm;
  
  if (field_time_bins_ == Teuchos::null){
    auto cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetGhosted(false);
    cvs->SetMesh(result->Map().Mesh());
    for (auto it=result->begin(); it!=result->end(); ++it){
      cvs->SetComponent(*it, result->Location(*it), number_bins_);
    }

    field_time_bins_ = Teuchos::rcp(new CompositeVector(*cvs));

    current_bin_ = 0;
    bin_time_[current_bin_] = 0.;
    update_time_ = S->time();        
  }
  
  const double dt = *S -> GetScalarData("dt");
  if (bin_time_[current_bin_] >= bin_period_){
    current_bin_ = (current_bin_ + 1)%number_bins_;
    for (auto it=result->begin(); it!=result->end(); ++it){
      Epetra_MultiVector& bin = *field_time_bins_ -> ViewComponent(*it, false);
      for (int j=0; j<bin.MyLength(); j++) bin[current_bin_][j] = 0.;
    }
    bin_time_[current_bin_] = 0.;      
  }
  

  bin_time_[current_bin_] += dt;
  for (auto it=result->begin(); it!=result->end(); ++it){
    Epetra_MultiVector& bin = *field_time_bins_ -> ViewComponent(*it, false);
    const Epetra_MultiVector& data = *base_field -> ViewComponent(*it, false);

    for (int j=0; j<bin.MyLength(); j++) {
      bin[current_bin_][j] = bin[current_bin_][j] + data[0][j]*dt;
    }
  }


  
  result -> PutScalar(0.);
  double total_time=0.;
  
  for (auto it=result->begin(); it!=result->end(); ++it){
    Epetra_MultiVector& bin = *field_time_bins_ -> ViewComponent(*it, false);
    Epetra_MultiVector& res = *result -> ViewComponent(*it, false);

  field_time_bins_ -> ViewComponent(*it, false)->Norm2(norm);
    
    for (int k=0; k<number_bins_; k++){
      //*vo_->os() << "k="<<k<<" bin_time_[k] "<<bin_time_[k]<<" norm "<< norm[k]<<"\n";      
      total_time += bin_time_[k];
      for (int j=0; j<bin.MyLength(); j++) {
        res[0][j] = res[0][j] + bin[k][j];
      }
    }

    if (total_time > 1e-12) {
      res.Scale(1./total_time);
    }
  }
}
  

bool TimeAverageFieldEvaluator::IsDependency(const Teuchos::Ptr<State>& S,
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

bool TimeAverageFieldEvaluator::ProvidesKey(Key key) const {
  return (key == my_key_);
}


void TimeAverageFieldEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
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



std::string
TimeAverageFieldEvaluator::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: time average" << std::endl;
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep != dependencies_.end(); ++dep) {
    result << "  Dep: " << *dep << std::endl;
  }
  result << std::endl;
  return result.str();
}

} // namespace
