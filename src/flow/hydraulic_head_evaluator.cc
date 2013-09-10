#include "hydraulic_head_evaluator.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziFlow {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
HydraulicHeadEvaluator::HydraulicHeadEvaluator(Teuchos::ParameterList& plist) :
    FieldEvaluator(plist) {

  my_key_ = plist.get<std::string>("evaluator name", "hydraulic_head");
  //my_mesh_ = "domain";

  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else if (my_key_ == std::string("hydraulic_head")) {
    my_mesh_ = "domain";
  } else if (my_key_.length() > std::string("hydraulic_head").length()) {
    my_mesh_ = my_key_.substr(0,my_key_.length() - std::string("hydraulic_head").length()-1);
  } else {
    ASSERT(0);
  }

  communicate_ = plist_.get<bool>("manage communication", false);
}

HydraulicHeadEvaluator::HydraulicHeadEvaluator(const HydraulicHeadEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_),
    my_mesh_(other.my_mesh_),
    communicate_(other.communicate_) {}


Teuchos::RCP<FieldEvaluator>
HydraulicHeadEvaluator::Clone() const {
  return Teuchos::rcp(new HydraulicHeadEvaluator(*this));
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void HydraulicHeadEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Require the field
  Teuchos::RCP<CompositeVectorFactory> my_fac = S->RequireField(my_key_, my_key_);

  if (!my_fac->owned()) {
    // requirements not yet set, set valid component
    S->RequireField(my_key_, "state")->SetMesh(S->GetMesh(my_mesh_))
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    // check plist for vis or checkpointing control
    bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  }
}


// ---------------------------------------------------------------------------
// Answers the question, has this Field changed since it was last requested
// for Field Key reqest.  Updates the field if needed.
// ---------------------------------------------------------------------------
bool HydraulicHeadEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) {
  UpdateField_(S);
  return true;
}


// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
inline
bool HydraulicHeadEvaluator::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {
  return false;  // no derivatives, though this should never be called
}


inline
bool HydraulicHeadEvaluator::IsDependency(const Teuchos::Ptr<State>& S, Key key) const {
  return false;
}

inline
bool HydraulicHeadEvaluator::ProvidesKey(Key key) const { return key == my_key_; }


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void HydraulicHeadEvaluator::UpdateField_(const Teuchos::Ptr<State>& S) {

  // NOTE: HydraulicHeadEvaluator data lives in state.
  Teuchos::RCP<CompositeVector> hh = S->GetFieldData(my_key_, "state");
  Epetra_MultiVector& hh_vec = *hh->ViewComponent("cell",false);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(my_mesh_);

  Teuchos::RCP<const CompositeVector> pressure = S->GetFieldData("pressure");
  const Epetra_MultiVector& pressure_vec = *pressure->ViewComponent("cell",false);
  
  hh_vec[0] = pressure_vec[0];

  // communicate if requested
  if (pressure->ghosted() && communicate_) {
    pressure->ScatterMasterToGhosted();
  }
}


} // namespace
} // namespace
