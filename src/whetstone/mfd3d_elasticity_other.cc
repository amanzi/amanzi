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

#include "mfd3d_elasticity.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* The stable discretization for Stokes: vectors at nodes plus normal
* components at faces. Fixed normal is used for the latter.
****************************************************************** */
int MFD3D_Elasticity::H1consistencyBernardiRaugel(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  int nrows = N.NumRows();
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int d = mesh_->space_dimension();
  AmanziGeometry::Point xv(d), tau(d), v1(d);

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xm = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // convolute tensors for non-zero modes
  std::vector<Tensor> vT, vKT;

  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      Tensor T(d, 2);
      T(i, j) = T(j, i) += 1.0;
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
  // 2D algorithm is separated out since it is fast
  N.PutScalar(0.0);

  if (d == 2) { 
    for (int n = 0; n < nnodes; n++) {
      int m = (n + 1) % nnodes;

      mesh_->node_get_coordinates(nodes[n], &xv);
      mesh_->node_get_coordinates(nodes[m], &v1);

      tau = v1 - xv;
      double length = norm(tau);
      tau /= length;

      AmanziGeometry::Point normal(d);
      normal[0] = tau[1];
      normal[1] =-tau[0];

      for (int i = 0; i < modes; i++) {
        v1 = vKT[i] * normal;
        double t = (tau * v1) * length / 2;

        N(2*n, i) += tau[0] * t;  
        N(2*n + 1, i) += tau[1] * t;

        N(2*m, i) += tau[0] * t;
        N(2*m + 1, i) += tau[1] * t;
      }
    }
  }

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);

    for (int i = 0; i < modes; i++) {
      v1 = vKT[i] * normal;
      double p = (normal * v1) / area;
      N(d*nnodes + n, i) += p * dirs[n];
    }
  }

  // calculate R coefM^{-1} R^T 
  DenseVector a1(modes), a2(modes), a3(modes);
  coefM.Inverse();

  for (int i = 0; i < nrows; i++) { 
    a1(0) = N(i, 0);  
    a1(1) = N(i, 1);  
    a1(2) = N(i, 2); 
    coefM.Multiply(a1, a3, false);

    for (int j = i; j < nrows; j++) {
      a2(0) = N(j, 0);  
      a2(1) = N(j, 1);  
      a2(2) = N(j, 2); 

      Ac(i, j) = a2 * a3;
    }
  }

  // calculate N (common algorihtm for 2D and 3D)
  N.PutScalar(0.0);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_coordinates(v, &xv);
    xv -= xm;

    // null modes
    int col(0);
    for (int k = 0; k < d; ++k) {
      N(d*n + k, col) = 1.0;
      col++;

      for (int l = k + 1; l < d; ++l) {
        N(d*n + k, col) = xv[l];  
        N(d*n + l, col) =-xv[k];
        col++;
      }
    }

    // non-null modes  
    for (int k = 0; k < d; ++k) {
      N(d*n + k, col) = xv[k];  
      col++;

      for (int l = k + 1; l < d; ++l) {
        N(d*n + k, col) = xv[l];
        N(d*n + l, col) = xv[k];
        col++;
      }
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    double area = mesh_->face_area(f);
    AmanziGeometry::Point normal(mesh_->face_normal(f));
    normal /= area;

    v1 = xf - xm;
    int m = d * nnodes + i;

    // null modes
    int col(0);
    for (int k = 0; k < d; ++k) {
      N(m, col) = normal[k];
      col++;

      for (int l = k + 1; l < d; ++l) {
        N(m, col) = v1[l] * normal[k] - v1[k] * normal[l];
        col++;
      }
    }

    // non-null modes
    for (int k = 0; k < d; ++k) {
      N(m, col) = v1[k] * normal[k];
      col++;

      for (int l = k + 1; l < d; ++l) {
        N(m, col) = v1[l] * normal[k] + v1[k] * normal[l];
        col++;
      }
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrixBernardiRaugel(
    int c, const Tensor& K, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nrows = A.NumRows();

  DenseMatrix N(nrows, d * (d + 1));

  int ok = H1consistencyBernardiRaugel(c, K, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Divergence matrix: vectors at nodes, normal components on faces.
* Fixed normal vector is used for the latter.
****************************************************************** */
int MFD3D_Elasticity::DivergenceMatrixBernardiRaugel(int c, DenseMatrix& A)
{
  int d = mesh_->space_dimension();

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  // mesh_->cell_get_nodes(c, &nodes);
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  // int n1 = d * nodes.size();
  // for (int n = 0; n < n1; ++n) A(0, n) = 0.0;

  for (int n = 0; n < nfaces; ++n) {
    double area = mesh_->face_area(faces[n]);
    // A(0, n1 + n) = area * dirs[n]; 
    A(0, n) = area * dirs[n]; 
  } 

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



