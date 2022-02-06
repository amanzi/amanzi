/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator for an changing cell volume.
*/

#include "errors.hh"
#include "Mesh.hh"

#include "EvaluatorDeformingCellVolume.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorDeformingCellVolume::EvaluatorDeformingCellVolume(Teuchos::ParameterList& plist)
    : EvaluatorSecondary(plist) {
  my_key_ = plist.get<std::string>("evaluator name", "cell_volume");

  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else if (my_key_ == std::string("cell_volume")) {
    my_mesh_ = "domain";
  } else if (my_key_.length() > std::string("cell_volume").length()) {
    my_mesh_ = my_key_.substr(0, my_key_.length() - std::string("cell_volume").length() - 1);
  } else {
    AMANZI_ASSERT(0);
  }

  // stick in the deformation key as my leaf node.
  Key deformation_key = plist.get<std::string>("deformation key");
  dependencies_.insert(std::make_pair(deformation_key, Tags::DEFAULT));
}


// ---------------------------------------------------------------------------
// Virutal copy constructor for Evaluator.
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator> EvaluatorDeformingCellVolume::Clone() const {
  return Teuchos::rcp(new EvaluatorDeformingCellVolume(*this));
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void EvaluatorDeformingCellVolume::EnsureCompatibility(State& S) {
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key_, Tags::DEFAULT);

  if (!my_fac.Owned()) {
    // requirements not yet set, claim ownership and set valid component
    S.Require<CompositeVector, CompositeVectorSpace>(my_key_, Tags::DEFAULT)
        .SetMesh(S.GetMesh(my_mesh_))->SetComponent("cell", AmanziMesh::CELL, 1);

    // check plist for vis or checkpointing control
    bool io_my_key =
        plist_.get<bool>(std::string("visualize ") + my_key_, true);
    std::cout << "vis " << my_key_ << "? " << io_my_key << std::endl;
    S.GetRecordW(my_key_, my_key_).set_io_vis(io_my_key);
    bool flag = plist_.get<bool>(std::string("checkpoint ") + my_key_, false);
    S.GetRecordW(my_key_, my_key_).set_io_vis(flag);
  }
}


// ---------------------------------------------------------------------------
// Evaluates the cell volume from the mesh values.
// ---------------------------------------------------------------------------
void EvaluatorDeformingCellVolume::Update_(State& S) {
  auto& cv = *S.GetW<CompositeVector>(my_key_, Tags::DEFAULT, my_key_).ViewComponent("cell");

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S.GetMesh(my_mesh_);

  // initialize from mesh
  int ncells = cv.MyLength();
  for (int c = 0; c != ncells; ++c) {
    cv[0][c] = mesh->cell_volume(c);
    if (cv[0][c] < 0.0)
      std::cout << "NEGATIVE CELL VOLUME cell " << c << ": " << cv[0][c] << std::endl;
  }
}

} // namespace Amanzi
