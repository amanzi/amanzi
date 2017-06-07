/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method on 2D surfaces in 3D.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of
* fluxes for a non-flat surface. Only the upper triangular part of 
* Wc is calculated. Darcy flux is scaled by the area!
* WARNING: routine works for scalar T only. 
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverseSurface(
    int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.NumRows()) return num_faces;  // matrix was not reshaped

  int dir, d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(c);

  // calculate cell normal
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf1 = mesh_->face_centroid(faces[0]);
  const AmanziGeometry::Point& xf2 = mesh_->face_centroid(faces[1]);
  AmanziGeometry::Point v1(d), v2(d), v3(d);

  v1 = (xf1 - xc)^(xf2 - xc);
  v1 /= norm(v1);

  // calculate projector
  Tensor P(d, 2); 
  for (int i = 0; i < d; i++) {
    P(i, i) = 1.0;
    for (int j = 0; j < d; j++) { 
      P(i, j) -= v1[i] * v1[j];
    }
  }

  // cell-based coordinate system
  v2 = xf1 - xc;
  v2 /= norm(v2);
  v3 = v1^v2; 

  // define new tensor
  Tensor PTP(d, 2);
  PTP = P * T * P;

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_face_normal(f, c);

    v1 = PTP * normal;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_face_normal(f, c);
      Wc(i, j) = (v1 * v2) / volume;
    }
  }

  // calculate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);

    R(i, 0) = v2 * (fm - cm);
    R(i, 1) = v3 * (fm - cm);
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy inverse mass matrix for surface: the standard algorithm
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseSurface(
    int cell, const Tensor& permeability, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d - 1);

  int ok = L2consistencyInverseSurface(cell, permeability, R, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, R, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Exterior normal to 2D face in 3D space.
****************************************************************** */
AmanziGeometry::Point MFD3D_Diffusion::mesh_face_normal(int f, int c)
{
  std::vector<AmanziGeometry::Point> vs;
  mesh_->face_get_coordinates(f, &vs);

  AmanziGeometry::Point tau(vs[1] - vs[0]);
  AmanziGeometry::Point normal = vs[0] - mesh_->cell_centroid(c);

  // orthogonalize and rescale normal
  double len = norm(tau);
  double s = (normal * tau) / len / len;
  normal -= s * tau;
  normal *= len / norm(normal);

  return normal;
}

}  // namespace WhetStone
}  // namespace Amanzi



