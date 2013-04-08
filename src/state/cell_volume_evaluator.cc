/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging cell volume.

------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "cell_volume_evaluator.hh"
#include "Mesh.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,CellVolumeEvaluator>
CellVolumeEvaluator::fac_("cell volume");

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
CellVolumeEvaluator::CellVolumeEvaluator(Teuchos::ParameterList& plist) :
    FieldEvaluator(plist),
    computed_once_(false) {

  my_key_ = plist.get<std::string>("evaluator name", "cell_volume");
  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else if (my_key_ == std::string("cell_volume")) {
    my_mesh_ = "domain";
  } else if (my_key_.length() > std::string("cell_volume").length()) {
    my_mesh_ = my_key_.substr(0,my_key_.length() - std::string("cell_volume").length()-1);
  } else {
    ASSERT(0);
  }

  communicate_ = plist_.get<bool>("manage communication", false);

}

CellVolumeEvaluator::CellVolumeEvaluator(const CellVolumeEvaluator& other) :
    FieldEvaluator(other),
    my_key_(other.my_key_),
    my_mesh_(other.my_mesh_),
    computed_once_(false),
    communicate_(other.communicate_) {}


Teuchos::RCP<FieldEvaluator>
CellVolumeEvaluator::Clone() const {
  return Teuchos::rcp(new CellVolumeEvaluator(*this));
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void CellVolumeEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Require the field
  Teuchos::RCP<CompositeVectorFactory> my_fac = S->RequireField(my_key_);

  if (!my_fac->owned()) {
    // requirements not yet set, claim ownership and set valid component
    S->RequireField(my_key_, my_key_)->SetMesh(S->GetMesh(my_mesh_))
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
bool CellVolumeEvaluator::HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) {
  if (!computed_once_) {
    // field DOES have to be computed at least once, even if it never changes.
    UpdateField_(S);
    computed_once_ = true;
    return true;
  }

  return false;
}


// ---------------------------------------------------------------------------
// Answers the question, Has This Field's derivative with respect to Key
// wrt_key changed since it was last requested for Field Key reqest.
// Updates the derivative if needed.
// ---------------------------------------------------------------------------
inline
bool CellVolumeEvaluator::HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
        Key request, Key wrt_key) {
  return false;  // no derivatives, though this should never be called
}


inline
bool CellVolumeEvaluator::IsDependency(const Teuchos::Ptr<State>& S, Key key) const {
  return false;
}

inline
bool CellVolumeEvaluator::ProvidesKey(Key key) const { return key == my_key_; }


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void CellVolumeEvaluator::UpdateField_(const Teuchos::Ptr<State>& S) {
  // NOTE: CellVolumeEvaluator owns its own data.
  Teuchos::RCP<CompositeVector> cv = S->GetFieldData(my_key_, my_key_);
  Epetra_MultiVector& cv_vec = *cv->ViewComponent("cell",false);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(my_mesh_);

  // initialize from mesh
  int ncells = cv_vec.MyLength();
  for (int c=0; c!=ncells; ++c) {
    cv_vec[0][c] = mesh->cell_volume(c);
  }

  // communicate if requested
  if (cv->ghosted() && communicate_) {
    cv->ScatterMasterToGhosted();
  }
}


} // namespace
