/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an changing cell volume.

------------------------------------------------------------------------- */
#include "errors.hh"
#include "Mesh.hh"

#include "deforming_cell_volume_evaluator.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
DeformingCellVolumeEvaluator::DeformingCellVolumeEvaluator(
        Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  my_key_ = plist.get<std::string>("evaluator name", "cell_volume");

  my_mesh_ = Keys::getDomain(my_key_);
  my_mesh_ = plist.get<std::string>("mesh name", my_mesh_);
  if (my_mesh_.empty()) {
    my_mesh_ = "domain";
  }

  // stick in the deformation key as my leaf node.
  Key deformation_key = Keys::readKey(plist, Keys::getDomain(my_key_), "deformation key", "base_porosity");
  dependencies_.insert(deformation_key);
}


// ---------------------------------------------------------------------------
// Copy constructor.
// ---------------------------------------------------------------------------
DeformingCellVolumeEvaluator::DeformingCellVolumeEvaluator(
        const DeformingCellVolumeEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    my_mesh_(other.my_mesh_) {}


// ---------------------------------------------------------------------------
// Virutal copy constructor for FieldEvaluator.
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator>
DeformingCellVolumeEvaluator::Clone() const {
  return Teuchos::rcp(new DeformingCellVolumeEvaluator(*this));
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void DeformingCellVolumeEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // Require the field
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_);
 
  if (!my_fac->Owned()) {
    // requirements not yet set, claim ownership and set valid component
    S->RequireField(my_key_, my_key_)->SetMesh(S->GetMesh(my_mesh_))
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    // check plist for vis or checkpointing control
    bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
    std::cout << "vis " << my_key_ << "? " << io_my_key << std::endl;
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  }
}


// ---------------------------------------------------------------------------
// Evaluates the cell volume from the mesh values.
// ---------------------------------------------------------------------------
void DeformingCellVolumeEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  // NOTE: CellVolumeEvaluator owns its own data.
  Epetra_MultiVector& cv = *S->GetFieldData(my_key_, my_key_)
      ->ViewComponent("cell",false);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(my_mesh_);

  // initialize from mesh
  int ncells = cv.MyLength();
  for (int c=0; c!=ncells; ++c) {
    cv[0][c] = mesh->cell_volume(c);
    if (cv[0][c] < 0.)
      std::cout << "NEGATIVE CELL VOLUME cell " << c << ": " << cv[0][c] << std::endl;
  }
}


void DeformingCellVolumeEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  result->PutScalar(0.);
  // Errors::Message message("Deforming cell volume's derivatives are not implemented");
  // Exceptions::amanzi_throw(message);
}


} // namespace
