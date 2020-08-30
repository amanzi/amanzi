/*
  WhetStone, Version 2.2
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

  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  R.Reshape(nfaces, d_ - 1);
  Wc.Reshape(nfaces, nfaces);

  double volume = mesh_->cell_volume(c);

  // calculate cell normal
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf1 = mesh_->face_centroid(faces[0]);
  const AmanziGeometry::Point& xf2 = mesh_->face_centroid(faces[1]);
  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);

  v1 = (xf1 - xc)^(xf2 - xc);
  v1 /= norm(v1);

  // calculate projector
  Tensor P(d_, 2); 
  for (int i = 0; i < d_; i++) {
    P(i, i) = 1.0;
    for (int j = 0; j < d_; j++) { 
      P(i, j) -= v1[i] * v1[j];
    }
  }

  // cell-based coordinate system
  v2 = xf1 - xc;
  v2 /= norm(v2);
  v3 = v1^v2; 

  // define new tensor
  Tensor PTP(d_, 2);
  PTP = P * T * P;

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_face_normal(f, c);

    v1 = PTP * normal;

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v4 = mesh_face_normal(f, c);
      Wc(i, j) = (v1 * v4) / volume;
    }
  }

  // calculate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);

    R(i, 0) = v2 * (fm - cm);
    R(i, 1) = v3 * (fm - cm);
  }

  return 0;
}


/* ******************************************************************
* Darcy inverse mass matrix for surface: the standard algorithm
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseSurface(
    int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverseSurface(c, K, R, W);
  if (ok) return ok;

  StabilityOptimized_(K, R, W);
  return 0;
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



