/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Release name: ara-to.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"

#include "mfd3d_diffusion.hh"
#include "tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of
* fluxes for a non-flat surface. Only the upper triangular part of 
* Wc is calculated. Darcy flux is scaled by the area!
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverseSurface(
    int cell, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.NumRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(cell);

  // calculate cell normal
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(cell);
  const AmanziGeometry::Point& xf1 = mesh_->face_centroid(faces[0]);
  const AmanziGeometry::Point& xf2 = mesh_->face_centroid(faces[1]);
  AmanziGeometry::Point cell_normal(d);
  cell_normal = (xf1 - xc)^(xf2 - xc);

  // calculate projector
  DenseMatrix P(d, d); 
  double a = L22(cell_normal);
  for (int i = 0; i < d; i++) {};

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = T * normal;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume);
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    for (int k = 0; k < d; k++) R(i, k) = fm[k] - cm[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



