/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d_diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* The stable discretization for Stokes: vectors at nodes plus normal
* components at faces.
****************************************************************** */
int MFD3D_Diffusion::H1consistencyNodeFace(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  int nrows = N.NumRows();
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int d = mesh_->space_dimension();
  AmanziGeometry::Point xv(d), tau(d), v1(d);

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nnodes, nfaces = faces.size();

  const AmanziGeometry::Point& xm = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // convolute tensors for non-zero modes
  std::vector<Tensor> vT, vKT;

  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      Tensor T(d, 2);
      T(i, j) = T(j, i) = 1.0;
      vT.push_back(T);

      Tensor KT(K * T);
      vKT.push_back(KT);
    }
  }

  // calculate exact integration matrix
  int modes = d * (d + 1) / 2;
  DenseMatrix coefM(modes, modes);

  for (int i = 0; i < modes; ++i) {
    for (int j = i; j < modes; ++j) {
      coefM(i, j) = DotTensor(vT[i], vKT[j]) * volume;
      coefM(j, i) = coefM(i, j);
    }
  }

  // to calculate matrix R, we use temporary matrix N 
  N.PutScalar(0.0);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    AmanziGeometry::Point normal(mesh_->face_normal(f));
    double area = mesh_->face_area(f);
    normal /= area;

    mesh_->face_get_nodes(f, &nodes);
    nnodes = nodes.size();

    if (d == 2) { 
      int m = (n + 1) % nnodes;
      tau[0] = normal[1];
      tau[1] =-normal[0];

      for (int i = 0; i < modes; i++) {
        v1 = vKT[i] * normal;
        double t = (tau * v1) * area / 2;
        double p = (normal * v1) * area;

        N(2*n, i) += tau[0] * t;  
        N(2*n + 1, i) += tau[1] * t;
        N(2*nnodes + n, i) += p;

        N(2*m, i) += tau[0] * t;
        N(2*m + 1, i) += tau[1] * t;
      }
    }
  }

  // calculate R K R^T / volume
  DenseVector a1(modes), a2(modes), a3(modes);
  coefM.Inverse();

  for (int i = 0; i < nrows; i++) { 
    a1(0) = N(i, 0);  
    a1(1) = N(i, 1);  
    a1(2) = N(i, 2); 

    for (int j = i; j < nrows; j++) {
      a2(0) = N(j, 0);  
      a2(1) = N(j, 1);  
      a2(2) = N(j, 2); 

      coefM.Multiply(a2, a3, false);
      Ac(i, j) = a1 * a3;
    }
  }

  // calculate N (common algorihtm for 2D and 3D)
  N.PutScalar(0.0);

  mesh_->cell_get_nodes(c, &nodes);
  nnodes = nodes.size();

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_coordinates(v, &xv);
    xv -= xm;

    // null mods
    int col(0);
    for (int k = 0; k < d; ++k) {
      N(2*n + k, col) = 1.0;
      col++;

      for (int l = k + 1; l < d; ++l) {
        N(2*n + k, col) = xv[l];  
        N(2*n + l, col) =-xv[k];
        col++;
      }
    }

    // non-null modes  
    for (int k = 0; k < d; ++k) {
      N(2*n + k, col) = xv[k];  
      col++;

      for (int l = k + 1; l < d; ++l) {
        N(2*n + k, col) = xv[l];
        N(2*n + k, col) = xv[k];
        col++;
      }
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    v1 = xf - xm;
    int l = d * nnodes + i;

    /*
    a1(l) = normal.x;
    a2(l) = p2D(xe.y, -xe.x) * normal;
    a3(l) = normal.y;

    a4(l) = xe.x * normal.x;
    a5(l) = p2D(xe.y, xe.x) * normal;
    a6(l) = xe.y * normal.y;
    */
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



