/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <cmath>
#include <set>
#include <vector>

#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "errors.hh"
#include "lapack.hh"

#include "ReconstructionCellCochain.hh"
#include "ReconstructionUtils.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructor allocates d + d^2 components in cells.
****************************************************************** */
ReconstructionCellCochain::ReconstructionCellCochain(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
  : d_(mesh->getSpaceDimension()), mesh_(mesh)
{
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_)->SetGhosted(true)->SetComponent(
    "cell", AmanziMesh::Entity_kind::CELL, d_ + d_ * d_);
  data_ = Teuchos::rcp(new CompositeVector(cvs, true));
  data_c_ = data_->ViewComponent("cell", true);
}


/* ******************************************************************
* Equality constrained least squares solved via the KKT system:
*   [ A^T A  C^T ] [u]   [A^T b]
*   [ C      0   ] [l] = [d    ]
* Coordinates are scaled by the cell size h and rows by the face area.
****************************************************************** */
void
ReconstructionCellCochain::Compute(const CompositeVector& flux)
{
  flux.ScatterMasterToGhosted("face");
  const auto& flux_f = *flux.ViewComponent("face", true);
  auto& data = *data_->ViewComponent("cell", true);

  int ncells_owned = mesh_->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nu = d_ + d_ * d_;

  std::set<AmanziMesh::Entity_ID> cells;
  WhetStone::DenseVector row(nu);

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& xc = mesh_->getCellCentroid(c);
    double h = std::pow(mesh_->getCellVolume(c), 1.0 / d_);

    const auto& own_faces = mesh_->getCellFaces(c);
    int nown = own_faces.size();
    int nk = nu + nown;

    getCellAllAdjCells(*mesh_, c, cells);
    cells.insert(c);

    std::set<AmanziMesh::Entity_ID> faces;
    for (int c1 : cells) {
      for (int f : mesh_->getCellFaces(c1)) faces.insert(f);
    }

    WhetStone::DenseMatrix K(nk, nk);
    WhetStone::DenseVector rhs(nk);
    K.PutScalar(0.0);
    rhs.PutScalar(0.0);

    for (int f : faces) {
      const auto& normal = mesh_->getFaceNormal(f);
      const auto& xf = mesh_->getFaceCentroid(f);
      double s = norm(normal);

      for (int i = 0; i < d_; ++i) {
        row(i) = normal[i] / s;
        for (int j = 0; j < d_; ++j) row(d_ + i * d_ + j) = row(i) * (xf[j] - xc[j]) / h;
      }
      double b = flux_f[0][f] / s;

      int k = 0;
      while (k < nown && own_faces[k] != f) ++k;

      if (k < nown) {
        for (int i = 0; i < nu; ++i) {
          K(nu + k, i) = K(i, nu + k) = row(i);
        }
        rhs(nu + k) = b;
      } else {
        for (int i = 0; i < nu; ++i) {
          for (int j = 0; j < nu; ++j) K(i, j) += row(i) * row(j);
          rhs(i) += row(i) * b;
        }
      }
    }

    int info, nrhs = 1;
    std::vector<int> ipiv(nk);
    WhetStone::DGESV_F77(&nk, &nrhs, K.Values(), &nk, ipiv.data(), rhs.Values(), &nk, &info);
    if (info != 0) {
      Errors::Message msg;
      msg << "Reconstruction patch of cell " << c << " (gid="
          << mesh_->getMap(AmanziMesh::CELL, false).GID(c)
          << ") does not determine a linear field; the mesh is too thin or degenerate.";
      Exceptions::amanzi_throw(msg);
    }

    for (int i = 0; i < d_; ++i) data[i][c] = rhs(i);
    for (int i = 0; i < d_ * d_; ++i) data[d_ + i][c] = rhs(d_ + i) / h;
  }

  data_->ScatterMasterToGhosted("cell");
}


/* ******************************************************************
* Reconstructed vector.
****************************************************************** */
AmanziGeometry::Point
ReconstructionCellCochain::getMean(int c) const
{
  AmanziGeometry::Point v(d_);
  for (int i = 0; i < d_; ++i) v[i] = (*data_c_)[i][c];
  return v;
}


AmanziGeometry::Point
ReconstructionCellCochain::getValue(int c, const AmanziGeometry::Point& p) const
{
  const auto& xc = mesh_->getCellCentroid(c);
  auto v = getMean(c);
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) v[i] += getGradient(c, i, j) * (p[j] - xc[j]);
  }
  return v;
}

} // namespace Operators
} // namespace Amanzi
