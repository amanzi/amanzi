/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  Miscalleneous collection of simple non-member functions.
*/

#include "PK_Utils.hh"

namespace Amanzi {

/* ******************************************************************
* Average permeability tensor in horizontal direction.
****************************************************************** */
void
PKUtils_CalculatePermeabilityFactorInWell(const Teuchos::Ptr<State>& S,
                                          Teuchos::RCP<Epetra_MultiVector>& Kxy)
{
  if (!S->HasRecord("permeability", Tags::DEFAULT)) return;

  const auto& cv = S->Get<CompositeVector>("permeability", Tags::DEFAULT);
  cv.ScatterMasterToGhosted("cell");
  const auto& perm = *cv.ViewComponent("cell", true);

  int ncells_wghost =
    S->GetMesh()->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  int dim = perm.NumVectors();

  Kxy = Teuchos::rcp(new Epetra_Vector(S->GetMesh()->getMap(AmanziMesh::Entity_kind::CELL, true)));

  for (int c = 0; c < ncells_wghost; c++) {
    (*Kxy)[0][c] = 0.0;
    int idim = std::max(1, dim - 1);
    for (int i = 0; i < idim; i++) (*Kxy)[0][c] += perm[i][c];
    (*Kxy)[0][c] /= idim;
  }
}


/* ******************************************************************
* Return coordinate of mesh entity (
****************************************************************** */
AmanziGeometry::Point
PKUtils_EntityCoordinates(int id, AmanziMesh::Entity_ID kind, const AmanziMesh::Mesh& mesh)
{
  if (kind == AmanziMesh::Entity_kind::FACE) {
    return mesh.getFaceCentroid(id);
  } else if (kind == AmanziMesh::Entity_kind::CELL) {
    return mesh.getCellCentroid(id);
  } else if (kind == AmanziMesh::Entity_kind::NODE) {
    int d = mesh.getSpaceDimension();
    AmanziGeometry::Point xn(d);
    xn = mesh.getNodeCoordinate(id);
    return xn;
  } else if (kind == AmanziMesh::Entity_kind::EDGE) {
    return mesh.getEdgeCentroid(id);
  }
  return AmanziGeometry::Point();
}


/* ******************************************************************
* Flux to full gradient
****************************************************************** */
void
PKUtils_FluxToVector(const State& S, const CompositeVector& flux, CompositeVector& grad)
{
  auto mesh = S.GetMesh();
  int d = mesh->getSpaceDimension();

  flux.ScatterMasterToGhosted("face");

  const auto& flux_f = *flux.ViewComponent("face", true);
  auto& grad_c = *grad.ViewComponent("cell");

  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    double vol = mesh->getCellVolume(c);
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    const auto& [faces, dirs] = mesh->getCellFacesAndDirections(c);
    int nfaces = faces.size();

    AmanziGeometry::Point tmp(d);
    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];
      const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
      const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
      tmp += normal * (flux_f[0][f] * dirs[i]);
    }
    tmp /= vol;

    for (int i = 0; i < d; ++i) grad_c[i][c] = tmp[i];
  }
}

} // namespace Amanzi
