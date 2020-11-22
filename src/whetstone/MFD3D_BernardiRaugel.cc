/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <tuple>
#include <vector>

#include "MeshLight.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_BernardiRaugel.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> MFD3D_BernardiRaugel::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::NODE, DOF_Type::POINT, d_));
  items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::NORMAL_COMPONENT, 1));
  return items;
}


/* ******************************************************************
* The stable discretization for Stokes: vectors at nodes plus normal
* components at faces. Fixed normal is used for the latter.
****************************************************************** */
int MFD3D_BernardiRaugel::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  AmanziGeometry::Point xv(d_), tau(d_), v1(d_);

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  int nrows = d_ * nnodes + nfaces;
  int nd = d_ * (d_ + 1);
  N.Reshape(nrows, nd);
  Ac.Reshape(nrows, nrows);

  const AmanziGeometry::Point& xm = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // convolute tensors for non-zero modes
  std::vector<Tensor> vT, vKT;

  for (int i = 0; i < d_; ++i) {
    for (int j = i; j < d_; ++j) {
      Tensor T(d_, 2);
      T(i, j) = T(j, i) = 1.0;
      vT.push_back(T);

      Tensor KT(K * T);
      vKT.push_back(KT);
    }
  }

  // calculate exact integration matrix
  int modes = d_ * (d_ + 1) / 2;
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

  if (d_ == 2) { 
    for (int n = 0; n < nnodes; n++) {
      int m = (n + 1) % nnodes;

      mesh_->node_get_coordinates(nodes[n], &xv);
      mesh_->node_get_coordinates(nodes[m], &v1);

      tau = v1 - xv;
      double length = norm(tau);
      tau /= length;

      AmanziGeometry::Point normal(d_);
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
      N(d_ * nnodes + n, i) += p * dirs[n];
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
    for (int k = 0; k < d_; ++k) {
      N(d_ * n + k, col) = 1.0;
      col++;

      for (int l = k + 1; l < d_; ++l) {
        N(d_ * n + k, col) = xv[l];  
        N(d_ * n + l, col) =-xv[k];
        col++;
      }
    }

    // non-null modes  
    for (int k = 0; k < d_; ++k) {
      N(d_ * n + k, col) = xv[k];  
      col++;

      for (int l = k + 1; l < d_; ++l) {
        N(d_ * n + k, col) = xv[l];
        N(d_ * n + l, col) = xv[k];
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
    int m = d_ * nnodes + i;

    // null modes
    int col(0);
    for (int k = 0; k < d_; ++k) {
      N(m, col) = normal[k];
      col++;

      for (int l = k + 1; l < d_; ++l) {
        N(m, col) = v1[l] * normal[k] - v1[k] * normal[l];
        col++;
      }
    }

    // non-null modes
    for (int k = 0; k < d_; ++k) {
      N(m, col) = v1[k] * normal[k];
      col++;

      for (int l = k + 1; l < d_; ++l) {
        N(m, col) = v1[l] * normal[k] + v1[k] * normal[l];
        col++;
      }
    }
  }

  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_BernardiRaugel::StiffnessMatrix(
    int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Advection matrix depends on velocity u.
****************************************************************** */
int MFD3D_BernardiRaugel::AdvectionMatrix(
    int c, const std::vector<AmanziGeometry::Point>& u, DenseMatrix& A)
{
  AMANZI_ASSERT(d_ == 2);

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  // calculate corner normals and weigths
  std::vector<double> w(nnodes, 0.0);
  AmanziGeometry::Point xv(d_);
  std::vector<AmanziGeometry::Point> N(nnodes, xv);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
   
    int j = (i + nfaces + 1) % nfaces;
    N[i] += normal * dirs[i] / 2.0;
    N[j] += normal * dirs[i] / 2.0;

    int v = nodes[i];
    mesh_->node_get_coordinates(v, &xv);
    double tmp = ((xv - xc) * normal) * dirs[i] / 6.0;
    w[i] += tmp;
    w[j] += tmp;
  }

  // populate matrix
  int ndofs = d_ * nnodes + nfaces;
  A.Reshape(ndofs, ndofs);
  A.PutScalar(0.0);
 
  for (int i = 0; i < nnodes; ++i) {
    int i1 = 2 * i;
    for (int j = 0; j < nnodes; ++j) {
      int j1 = 2 * j;
      A(i1, j1) = w[i] * u[i][0] * N[j][0];
      A(i1, j1 + 1) = w[i] * u[i][1] * N[j][0];
      A(i1 + 1, j1) = w[i] * u[i][0] * N[j][1];
      A(i1 + 1, j1 + 1) = w[i] * u[i][1] * N[j][1];
    }    
  } 

  return 0;
}


/* ******************************************************************
* Divergence matrix: vectors at nodes, normal components on faces.
* Fixed normal vector is used for the latter.
****************************************************************** */
int MFD3D_BernardiRaugel::DivergenceMatrix(int c, DenseMatrix& A)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  int n1 = d_ * nodes.size();
  A.Reshape(1, n1 + nfaces);

  for (int n = 0; n < n1; ++n) A(0, n) = 0.0;

  for (int n = 0; n < nfaces; ++n) {
    double area = mesh_->face_area(faces[n]);
    A(0, n1 + n) = area * dirs[n]; 
  } 

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi



