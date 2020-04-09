/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "Mesh.hh"
#include "errors.hh"

#include "deforming_cell_volume_evaluator.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
DeformingCellVolumeEvaluator::DeformingCellVolumeEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype(plist)
{
  my_key_ = plist.get<std::string>("evaluator name", "cell_volume");

  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else if (my_key_ == std::string("cell_volume")) {
    my_mesh_ = "domain";
  } else if (my_key_.length() > std::string("cell_volume").length()) {
    my_mesh_ = my_key_.substr(
      0, my_key_.length() - std::string("cell_volume").length() - 1);
  } else {
    AMANZI_ASSERT(0);
  }

  // stick in the deformation key as my leaf node.
  Key deformation_key = plist.get<std::string>("deformation key");
  dependencies_.emplace_back(deformation_key);
}

// ---------------------------------------------------------------------------
// Copy constructor.
// ---------------------------------------------------------------------------
DeformingCellVolumeEvaluator::DeformingCellVolumeEvaluator(
  const DeformingCellVolumeEvaluator& other)
  : EvaluatorSecondaryMonotype(other), my_mesh_(other.my_mesh_)
{}

// ---------------------------------------------------------------------------
// Virutal copy constructor for Evaluator.
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
DeformingCellVolumeEvaluator::Clone() const
{
  return Teuchos::rcp(new DeformingCellVolumeEvaluator(*this));
}

// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void
DeformingCellVolumeEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // Require the field
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_);

  if (!my_fac->Owned()) {
    // requirements not yet set, claim ownership and set valid component
    S->RequireField(my_key_, my_key_)
      .SetMesh(S->GetMesh(my_mesh_))
      .SetComponent("cell", AmanziMesh::CELL, 1);

    // check plist for vis or checkpointing control
    bool io_my_key =
      plist_.get<bool>(std::string("visualize ") + my_key_, true);
    std::cout << "vis " << my_key_ << "? " << io_my_key << std::endl;
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key =
      plist_.get<bool>(std::string("checkpoint ") + my_key_, false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);
  }
}

// ---------------------------------------------------------------------------
// Evaluates the cell volume from the mesh values.
// ---------------------------------------------------------------------------
void
DeformingCellVolumeEvaluator::EvaluateField_(
  const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  // NOTE: CellVolumeEvaluator owns its own data.
  MultiVector_type& cv =
    *S->GetFieldData(my_key_, my_key_)->ViewComponent("cell", false);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh(my_mesh_);

  // initialize from mesh
  int ncells = cv.getLocalLength();
  for (int c = 0; c != ncells; ++c) {
    cv[0][c] = mesh->cell_volume(c, false);
    if (cv[0][c] < 0.)
      std::cout << "NEGATIVE CELL VOLUME cell " << c << ": " << cv[0][c]
                << std::endl;
  }
}

void
DeformingCellVolumeEvaluator::EvaluateFieldPartialDerivative_(
  const Teuchos::Ptr<State>& S, Key wrt_key,
  const Teuchos::Ptr<CompositeVector>& result)
{
  result->putScalar(0.);
  // Errors::Message message("Deforming cell volume's derivatives are not
  // implemented"); throw(message);
}

} // namespace Amanzi
