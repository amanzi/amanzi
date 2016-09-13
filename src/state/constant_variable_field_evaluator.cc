/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.



*/

#include "CompositeVectorFunctionFactory.hh"
#include "constant_variable_field_evaluator.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
ConstantVariableFieldEvaluator::ConstantVariableFieldEvaluator(Teuchos::ParameterList& plist) :
    FieldEvaluator(plist) {

  my_key_ = plist_.get<std::string>("evaluator name");

}


// ---------------------------------------------------------------------------
// Copy constructor
// ---------------------------------------------------------------------------
ConstantVariableFieldEvaluator::ConstantVariableFieldEvaluator(const ConstantVariableFieldEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_)
{}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator> ConstantVariableFieldEvaluator::Clone() const {
  return Teuchos::rcp(new ConstantVariableFieldEvaluator(*this));
}

// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
void ConstantVariableFieldEvaluator::operator=(const FieldEvaluator& other) {
#if ENABLE_DBC
  const ConstantVariableFieldEvaluator* other_p =
      dynamic_cast<const ConstantVariableFieldEvaluator*>(&other);
  ASSERT(other_p != NULL);
#else
  const ConstantVariableFieldEvaluator* other_p =
      static_cast<const ConstantVariableFieldEvaluator*>(&other);
#endif
  ASSERT(my_key_ == other_p->my_key_);

  requests_ = other_p->requests_;
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void ConstantVariableFieldEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Require the field and claim ownership.
  S->RequireField(my_key_, my_key_);
  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  
}

// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// ---------------------------------------------------------------------------
bool ConstantVariableFieldEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) {
  Teuchos::OSTab tab = vo_->getOSTab();

  return false;

}


// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
inline
bool ConstantVariableFieldEvaluator::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {

  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "INDEPENDENT Variable derivative requested by " << request
          << " has not changed." << std::endl;
  }
  return false;  // no derivatives, though this should never be called
}


inline
bool ConstantVariableFieldEvaluator::IsDependency(const Teuchos::Ptr<State>& S, Key key) const {
  return false;
}

inline
bool ConstantVariableFieldEvaluator::ProvidesKey(Key key) const { return key == my_key_; }


// ---------------------------------------------------------------------------
// String representation of this evaluator
// ---------------------------------------------------------------------------
std::string
ConstantVariableFieldEvaluator::WriteToString() const {
  std::stringstream result;
  result << my_key_ << std::endl
         << "  Type: independent" << std::endl
         << std::endl;
  return result.str();
}

}
