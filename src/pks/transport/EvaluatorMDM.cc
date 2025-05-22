/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#include "TensorVector.hh"
#include "MDMFactory.hh"
#include "MDMPartition.hh"
#include "EvaluatorMDM.hh"

namespace Amanzi {
namespace Transport {

EvaluatorMDM::EvaluatorMDM(Teuchos::ParameterList& plist) :
  EvaluatorSecondary(plist)
{
  std::string domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  is_surface_ = plist_.get<bool>("is surface", Keys::in(domain, "surface"));

  if (!is_surface_) {
    poro_key_ = Keys::readKey(plist_, domain, "porosity", "porosity");
    dependencies_.insert(KeyTag{poro_key_, tag});
  }

  if (is_surface_) {
    velocity_key_ = Keys::readKey(plist_, domain, "velocity", "water_velocity");
  } else {
    velocity_key_ = Keys::readKey(plist_, domain, "velocity", "darcy_velocity");
  }
  dependencies_.insert(KeyTag{velocity_key_, tag});
}

void
EvaluatorMDM::Update_(State& S)
{
  KeyTag key_tag = my_keys_.front();
  Tag tag = key_tag.second;

  TensorVector& D = S.GetW<TensorVector>(my_keys_.front().first, tag, my_keys_.front().first);
  const Epetra_MultiVector& velocity = *S.Get<CompositeVector>(velocity_key_, tag)
    .ViewComponent("cell", false);

  const AmanziMesh::Mesh& m = *S.GetMesh(Keys::getDomain(key_tag.first));

  const double& time = S.Get<double>("time", tag);
  int space_dim = m.getSpaceDimension();

  // NOTE: here we always assume the axi-symmetry is the vertical.  Amanzi
  // computes this by looking at permeability.
  int axi_symmetry = space_dim - 1;
  AmanziGeometry::Point v(space_dim);

  if (is_surface_) {
    for (int c = 0; c != D.size(); ++c) {
      for (int i = 0; i != space_dim; ++i) v[i] = velocity[i][c];
      D[c] = mdms_->second[(*mdms_->first)[c]]->mech_dispersion(
        time, m.getCellCentroid(c), v, axi_symmetry, 1.0, 1.0);
    }

  } else {
    const Epetra_MultiVector& poro = *S.Get<CompositeVector>(poro_key_, tag)
      .ViewComponent("cell", false);

    for (int c = 0; c != D.size(); ++c) {
      for (int i = 0; i != space_dim; ++i) v[i] = velocity[i][c];
      D[c] = mdms_->second[(*mdms_->first)[c]]->mech_dispersion(
        time, m.getCellCentroid(c), v, axi_symmetry, 1.0, poro[0][c]);
    }
  }
}


void
EvaluatorMDM::EnsureCompatibility(State& S)
{
  KeyTag my_keytag = my_keys_.front();
  auto& fac = S.Require<TensorVector, TensorVector_Factory>(my_keytag.first, my_keytag.second, my_keytag.first);
  EnsureCompatibility_Flags_(S);

  if (fac.map().Mesh() != Teuchos::null && mdms_ == Teuchos::null) {
    auto mesh = fac.map().Mesh();
    auto plist = Teuchos::rcpFromRef(plist_);
    bool flag(false);
    mdms_ = CreateMDMPartition(fac.map().Mesh(), Teuchos::sublist(plist, "mechanical dispersion parameters"), flag);

    if (!is_surface_) {
      S.Require<CompositeVector,CompositeVectorSpace>(poro_key_, my_keytag.second)
        .SetMesh(mesh)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }

    S.Require<CompositeVector,CompositeVectorSpace>(velocity_key_, my_keytag.second)
      .SetMesh(mesh)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, mesh->getSpaceDimension());

    CompositeVectorSpace space;
    space.SetMesh(mesh)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    fac.set_map(space);
  }

  for (const auto& dep : dependencies_) {
    S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
  }
}


} // namespace Transport
} // namespace Amanzi
