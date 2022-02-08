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
    : EvaluatorSecondary(plist)
{
  if (plist.isParameter("mesh name")) {
    my_mesh_ = plist.get<std::string>("mesh name");
  } else {
    my_mesh_ = Keys::getDomain(my_keys_.front().first);
  }

  // stick in the deformation key as my leaf node.
  Key deformation_key = plist.get<std::string>("deformation key");
  Tag deformation_tag(plist.get<std::string>("tag", ""));
  dependencies_.insert(KeyTag(deformation_key, deformation_tag));
}


// ---------------------------------------------------------------------------
// Virutal copy constructor for Evaluator.
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator> EvaluatorDeformingCellVolume::Clone() const {
  return Teuchos::rcp(new EvaluatorDeformingCellVolume(*this));
}


void EvaluatorDeformingCellVolume::EnsureCompatibility(State& S) {
  Key key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  S.Require<CompositeVector,CompositeVectorSpace>(key, tag)
    .SetMesh(S.GetMesh(my_mesh_))
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  EnsureCompatibility_Flags_(S);
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
