/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Tag.hh"
#include "State.hh"
#include "EvaluatorVelocityReconstruction.hh"

namespace Amanzi {

EvaluatorVelocityReconstruction::EvaluatorVelocityReconstruction(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace>(plist)
{
  std::string domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  flux_key_ = Keys::readKey(plist, domain, "water flux", "water_flux");
  dependencies_.insert({flux_key_, tag});

  molar_dens_key_ = Keys::readKey(plist, domain, "molar density liquid", "molar_density_liquid");
  dependencies_.insert({molar_dens_key_, tag});
}


Teuchos::RCP<Evaluator>
EvaluatorVelocityReconstruction::Clone() const {
  return Teuchos::rcp(new EvaluatorVelocityReconstruction(*this));
}


void
EvaluatorVelocityReconstruction::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& flux =
    *S.Get<CompositeVector>(flux_key_, tag).ViewComponent("face", true);
  const Epetra_MultiVector& nliq_c =
    *S.Get<CompositeVector>(molar_dens_key_, tag).ViewComponent("cell", false);

  Epetra_MultiVector& velocity =
    *results[0]->ViewComponent("cell", false);

  int d(mesh_->getSpaceDimension());
  AmanziGeometry::Point local_velocity(d);

  Teuchos::LAPACK<int, double> lapack;
  Teuchos::SerialDenseMatrix<int, double> matrix(d, d);
  double rhs[d];

  for (int c = 0; c != velocity.MyLength(); ++c) {
    auto faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int i = 0; i != d; ++i) rhs[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n = 0; n != nfaces; ++n) { // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
      double area = mesh_->getFaceArea(f);

      for (int i = 0; i != d; ++i) {
        rhs[i] += normal[i] * flux[0][f];
        matrix(i, i) += normal[i] * normal[i];
        for (int j = i + 1; j < d; ++j) { matrix(j, i) = matrix(i, j) += normal[i] * normal[j]; }
      }
    }

    int info;
    lapack.POSV('U', d, 1, matrix.values(), d, rhs, d, &info);

    for (int i = 0; i != d; ++i) velocity[i][c] = rhs[i] / nliq_c[0][c];
  }
}



void
EvaluatorVelocityReconstruction::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;
  CompositeVectorSpace& cvs = S.Require<CompositeVector,CompositeVectorSpace>(my_keys_.front().first, my_keys_.front().second);
  if (mesh_ == Teuchos::null && cvs.Mesh() != Teuchos::null) {
    mesh_ = cvs.Mesh();

    S.Require<CompositeVector,CompositeVectorSpace>(flux_key_, tag)
      .SetMesh(mesh_)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1)
      ->SetGhosted("true");

    S.Require<CompositeVector,CompositeVectorSpace>(molar_dens_key_, tag)
      .SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->SetGhosted("true");
  }
}

} // namespace Amanzi
