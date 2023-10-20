/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

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
  domain_ = Keys::readDomain(plist, "", Keys::getDomain(my_keys_.front().first));

  // stick in the deformation key as my leaf node.
  Key deformation_key = Keys::readKey(
    plist, Keys::getDomain(my_keys_.front().first), "deformation key", "base_porosity");
  Tag deformation_tag(plist.get<std::string>("dependency tag", my_keys_.front().second.get()));
  dependencies_.insert(KeyTag(deformation_key, deformation_tag));
}


// ---------------------------------------------------------------------------
// Virutal copy constructor for Evaluator.
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorDeformingCellVolume::Clone() const
{
  return Teuchos::rcp(new EvaluatorDeformingCellVolume(*this));
}


void
EvaluatorDeformingCellVolume::EnsureCompatibility(State& S)
{
  Key key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  S.Require<CompositeVector, CompositeVectorSpace>(key, tag, key)
    .SetMesh(S.GetMesh(domain_))
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // For dependencies, all we really care is whether there is an evaluator or
  // not.  We do not use the data at all.
  for (const auto& dep : dependencies_) { S.RequireEvaluator(dep.first, dep.second); }
  EnsureCompatibility_Flags_(S);
}


// ---------------------------------------------------------------------------
// Evaluates the cell volume from the mesh values.
// ---------------------------------------------------------------------------
void
EvaluatorDeformingCellVolume::Update_(State& S)
{
  auto key_tag = my_keys_.front();
  auto& cv =
    *S.GetW<CompositeVector>(key_tag.first, key_tag.second, key_tag.first).ViewComponent("cell");
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S.GetMesh(domain_);

  // initialize from mesh
  int ncells = cv.MyLength();
  for (int c = 0; c != ncells; ++c) {
    cv[0][c] = mesh->getCellVolume(c);
    if (cv[0][c] < 0.0)
      std::cout << "NEGATIVE CELL VOLUME cell " << c << ": " << cv[0][c] << std::endl;
  }
}

} // namespace Amanzi
