/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method.
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
* Consistency condition for inner product in space of fluxes. 
* Only upper triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
* Here R^T N = |c| K.
* Fluxes include face areas!
****************************************************************** */
int MFD3D_Diffusion::L2consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int nfaces = faces.size();
  if (nfaces != N.NumRows()) return nfaces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(c);

  AmanziGeometry::Point v1(d), v2(d);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  Tensor Kinv(K);
  Kinv.Inverse();
  Kinv.Transpose();

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    v2 = Kinv * (fm - cm);

    int i0 = (symmetry ? i : 0);
    for (int j = i0; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
      v1 = fm - cm;
      Mc(i, j) = (v1 * v2) / volume;
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    for (int k = 0; k < d; k++) N(i, k) = normal[k] * dirs[i];
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of the mass matrix in the space
* of Darcy fluxes. Only the upper triangular part of matrix 
* Wc = N (N^T R)^{-1} N^T is calculated. Here N^T R = |c| K.
* Flux is scaled by face area!
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverse(
    int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc, bool symmetry)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.NumRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(c);

  Tensor Kt(K);
  Kt.Transpose();

  // Since N is scaled by K, N = N0 * K, we us tensor K in the
  // inverse L2 consistency term.
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = Kt * normal;

    int i0 = (symmetry ? i : 0);
    for (int j = i0; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume);
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    for (int k = 0; k < d; k++) R(i, k) = fm[k] - cm[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction. 
* Only the upper triangular part of Ac is calculated.
* The degrees of freedom are at nodes.
****************************************************************** */
int MFD3D_Diffusion::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int num_nodes = nodes.size();
  if (num_nodes != N.NumRows()) return num_nodes;  // matrix was not reshaped

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(c);
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  /* to calculate matrix R, we use temporary matrix N */
  N.PutScalar(0.0);

  int num_faces = faces.size();
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d == 2) {
        u = 0.5 * dirs[i]; 
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        mesh_->node_get_coordinates(v, &p);
        mesh_->node_get_coordinates(vnext, &pnext);
        mesh_->node_get_coordinates(vprev, &pprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1^v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = FindPosition_(v, nodes);
      for (int k = 0; k < d; k++) N(pos, k) += normal[k] * u;
    }
  }

  for (int i = 0; i < num_nodes; i++) {  // calculate R K R^T / volume
    for (int k = 0; k < d; k++) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < num_nodes; j++) {
      for (int k = 0; k < d; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);
  for (int i = 0; i < num_nodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    for (int k = 0; k < d; k++) N(i, k) = p[k] - cm[k];
    N(i, d) = 1;  // additional column is added to the consistency condition
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction. 
* Only the upper triangular part of Ac is calculated.
* The degrees of freedom are on edges.
****************************************************************** */
int MFD3D_Diffusion::H1consistencyEdge(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List edges, faces, face_edges;
  std::vector<int> dirs, map;

  int d = mesh_->space_dimension();
  ASSERT(d == 2);

  mesh_->cell_get_edges(c, &edges);
  int num_edges = edges.size();
  if (num_edges != N.NumRows()) return num_edges;  // matrix was not reshaped

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int num_faces = faces.size();

  // to calculate matrix R, we use temporary matrix N 
  N.PutScalar(0.0);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    if (d == 2) {
      for (int k = 0; k < d; k++) N(i, k) += normal[k] * dirs[i];
    } else {
      mesh_->face_to_cell_edge_map(f, c, &map);
      int num_face_edges = map.size();

      for (int j = 0; j < num_face_edges; j++) {
        for (int k = 0; k < d; k++) N(map[j], k) += normal[k] * dirs[i];
      }
    }
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d), v2(d);
  double volume = mesh_->cell_volume(c);
  for (int i = 0; i < num_edges; i++) {
    for (int k = 0; k < d; k++) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < num_edges; j++) {
      for (int k = 0; k < d; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  // calculate N
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  for (int i = 0; i < num_edges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    for (int k = 0; k < d; k++) N(i, k) = xe[k] - xc[k];
    N(i, d) = 1;  // additional column is added to the consistency condition
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix in space of fluxes: the standard algorithm
****************************************************************** */
int MFD3D_Diffusion::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nfaces = M.NumRows();

  DenseMatrix N(nfaces, d);

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistency(c, Kinv, N, M, true);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix in space of fluxes for non-symmetric tensor
****************************************************************** */
int MFD3D_Diffusion::MassMatrixNonSymmetric(int c, const Tensor& K, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nfaces = M.NumRows();

  DenseMatrix N(nfaces, d);

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistency(c, Kinv, N, M, false);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalarNonSymmetric_(c, N, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse mass matrix in the space of fluxes: the standard algorithm
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, R, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse mass matrix for non-symmetric PD tensor.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseNonSymmetric(int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverse(c, K, R, W, false);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalarNonSymmetric_(c, R, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Diffusion::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nnodes = A.NumRows();

  DenseMatrix N(nnodes, d + 1);

  int ok = H1consistency(c, K, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Diffusion::StiffnessMatrixOptimized(int c, const Tensor& K, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nnodes = A.NumRows();

  DenseMatrix N(nnodes, d + 1);

  int ok = H1consistency(c, K, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(K, N, A);
  return ok;
}


/* ******************************************************************
* Stiffness matrix: the M-matrix approach
****************************************************************** */
int MFD3D_Diffusion::StiffnessMatrixMMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nnodes = A.NumRows();

  DenseMatrix N(nnodes, d + 1);

  int ok = H1consistency(c, K, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  // scaling of matrix A for numerical stability
  double s = A.Trace() / nnodes;
  A /= s;

  int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE;
  StabilityMMatrix_(c, N, A, objective);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  A *= s;
  simplex_functional_ *= s;
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. DOFs are on edges.
****************************************************************** */
int MFD3D_Diffusion::StiffnessMatrixEdge(int c, const Tensor& K, DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nedges = A.NumRows();

  DenseMatrix N(nedges, d + 1);

  int ok = H1consistencyEdge(c, K, N, A);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* *****************************************************************
* Recover gradient from solution, which is either the edge-based 
* fluxes of node-based pressures. The algorithm is common if both
* N and R are used. Here we use simplified versions.
***************************************************************** */
int MFD3D_Diffusion::RecoverGradient_MassMatrix(
    int c, const std::vector<double>& solution, AmanziGeometry::Point& gradient)
{
  int d = mesh_->space_dimension();
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int num_faces = faces.size();

  gradient.set(0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);

    for (int k = 0; k < d; k++) {
      double Rik = fm[k] - cm[k];
      gradient[k] += Rik * solution[i] * dirs[i];
    }
  }

  gradient *= -1.0 / mesh_->cell_volume(c);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* *****************************************************************
* Recover gradient from solution, which is either the edge-based 
* fluxes of node-based pressures. The algorithm is common if both
* N and R are used. Here we use simplified versions.
***************************************************************** */
int MFD3D_Diffusion::RecoverGradient_StiffnessMatrix(
    int c, const std::vector<double>& solution, AmanziGeometry::Point& gradient)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int num_nodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int num_faces = faces.size();

  // populate matrix R (should be a separate routine lipnikov@lanl.gv)
  int d = mesh_->space_dimension();
  DenseMatrix R(num_nodes, d);
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  R.PutScalar(0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d == 2) {
        u = 0.5 * dirs[i]; 
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        mesh_->node_get_coordinates(v, &p);
        mesh_->node_get_coordinates(vnext, &pnext);
        mesh_->node_get_coordinates(vprev, &pprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1^v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = FindPosition_(v, nodes);
      for (int k = 0; k < d; k++) R(pos, k) += normal[k] * u;
    }
  }

  gradient.set(0.0);
  for (int i = 0; i < num_nodes; i++) {
    for (int k = 0; k < d; k++) {
      gradient[k] += R(i, k) * solution[i];
    }
  }
  gradient *= 1.0 / mesh_->cell_volume(c);

  return 0;
}


/* *****************************************************************
*  OTHER ROUTINES
***************************************************************** */

/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of  
* fluxes. Only the upper triangular part of Wc is calculated.
* Flux is the integral average.
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverseScaled(
    int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.NumRows()) return num_faces;  // matrix was not reshaped

  // calculate areas of possibly curved faces
  std::vector<double> areas(num_faces, 0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->face_normal(f));
  }

  // populate matrix W_0
  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = K * normal;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    for (int k = 0; k < d; k++) R(i, k) = (fm[k] - cm[k]) * areas[i];
  }

  /* Internal verification 
  DenseMatrix NtR(d, d);
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      NtR(i, j) = 0.0;
      for (int k = 0; k < num_faces; k++) {
        const AmanziGeometry::Point& v1 = mesh_->face_normal(faces[k]);
        NtR(i, j) += v1[i] * R(k, j) / areas[k] * dirs[k];
      }
    }
  }
  */
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of 
* fluxes. Only the upper triangular part of Wc is calculated.
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverseDivKScaled(
    int c, const Tensor& K, double kmean, const AmanziGeometry::Point& kgrad,
    DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List faces, nodes;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.NumRows()) return num_faces;  // matrix was not reshaped

  // calculate areas of possibly curved faces
  std::vector<double> areas(num_faces, 0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->face_normal(f));
  }

  // populate matrix W_0
  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d), v2(d);
  double volume = mesh_->cell_volume(c);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = (K * normal) / kmean;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(c);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    if (d == 2) { 
      mesh_->face_get_nodes(f, &nodes);

      mesh_->node_get_coordinates(nodes[0], &v1);
      mesh_->node_get_coordinates(nodes[1], &v2);

      v1 -= cm;
      v2 -= cm;

      double k1 = kmean + kgrad * v1;
      double k2 = kmean + kgrad * v2;
      double km = k1 + k2;
      double tmp = areas[i] / 6;
      for (int k = 0; k < d; k++) {
        double vm = v1[k] + v2[k];
        R(i, k) = (k1 * v1[k] + k2 * v2[k] + km * vm) * tmp;
      }
    } else {
      ASSERT(false);
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse mass matrix in flux space via optimization, experimental.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseScaled(int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverseScaled(c, K, R, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
 
  StabilityScalar(c, R, W);
  RescaleMassMatrixInverse_(c, W);

  return ok;
}


/* ******************************************************************
* Mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseMMatrixHex(
    int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityMMatrixHex_(c, K, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix for a polyhedral element via simplex method.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseMMatrix(
    int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  // use boolean flag to populate the whole matrix
  int ok = L2consistencyInverse(c, K, R, W, false);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  // scaling of matrix W for numerical stability
  double s = W.Trace() / nfaces;
  W /= s;

  ok = StabilityMMatrix_(c, R, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  W *= s;
  simplex_functional_ *= s;
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Inverse mass matrix via optimization.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseOptimized(
    int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(K, R, W);
  return ok;
}


/* ******************************************************************
* Inverse mass matrix via optimization, experimental.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseOptimizedScaled(
    int c, const Tensor& K, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverseScaled(c, K, R, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
 
  ok = StabilityOptimized(K, R, W);
  RescaleMassMatrixInverse_(c, W);

  return ok;
}


/* ******************************************************************
* Inverse mass matrix via optimization, experimental.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseDivKScaled(
    int c, const Tensor& K, double kmean, const AmanziGeometry::Point& kgrad, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);

  int ok = L2consistencyInverseDivKScaled(c, K, kmean, kgrad, R, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
 
  StabilityScalar(c, R, W);
  RescaleMassMatrixInverse_(c, W);

  return ok;
}


/* ******************************************************************
* Rescale matrix to area-weighted fluxes.
****************************************************************** */
void MFD3D_Diffusion::RescaleMassMatrixInverse_(int c, DenseMatrix& W)
{
  Entity_ID_List faces;

  mesh_->cell_get_faces(c, &faces);
  int num_faces = faces.size();

  // calculate areas of possibly curved faces
  std::vector<double> areas(num_faces, 0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->face_normal(f));
  }

  // back to area-weighted fluxes
  for (int i = 0; i < num_faces; i++) {
    for (int j = 0; j < num_faces; j++) W(i, j) *= areas[i] * areas[j];
  }
}


/* ******************************************************************
* A simple monotone stability term for a 2D or 3D brick element. 
****************************************************************** */
int MFD3D_Diffusion::StabilityMMatrixHex_(int c, const Tensor& K, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = 2 * d;

  // symmetrize the consistency matrix
  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(j, i) = M(i, j);
  }

  // create groups of quasi-parallel faces
  int map[nrows];
  for (int i = 0; i < nrows; i++) map[i] = i;

  Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  int i1, i2, k, l;
  double s1, s2, area1, area2;
  for (int i = 0; i < d-1; i++) {
    int i1 = 2*i;
    int i2 = i1 + 1;

    int f = faces[i1];
    const AmanziGeometry::Point& normal1 = mesh_->face_normal(f);
    area1 = mesh_->face_area(f);

    s1 = 1.0;
    for (int j = i2; j < nrows; j++) {
      int f = faces[j];
      const AmanziGeometry::Point& normal2 = mesh_->face_normal(f);
      area2 = mesh_->face_area(f);

      s2 = (normal1 * normal2) * (dirs[i] * dirs[j]) / (area1 * area2);
      if (s2 < s1) {  // swap map values in positions i2 and j
        k = map[i2];
        map[i2] = j;
        map[j] = k;
        s1 = s2;
      }
    }
    // if (s1 >= 0.0) return -1;  // hex is too disturted
  }

  // define transformed tensor
  Tensor T1(d, 2);
  AmanziGeometry::Point areas(d);
  for (int i = 0; i < d; i++) {
    k = map[2*i];
    int f = faces[k];
    const AmanziGeometry::Point& normal1 = mesh_->face_normal(f);
    area1 = mesh_->face_area(f);
    areas[i] = area1;

    for (int j = i; j < d; j++) {
      l = map[2*j];
      f = faces[l];
      const AmanziGeometry::Point& normal2 = mesh_->face_normal(f);
      area2 = mesh_->face_area(f);

      s1 = (K * normal1) * normal2 * (dirs[k] * dirs[l]) / (area1 * area2);
      if (i-j) {
        T1(i, j) = T1(j, i) = -fabs(s1);
      } else {
        T1(i, i) = s1;
      }
    }
  }

  // verify SPD property
  if (d == 3) {
    double lower, upper;
    T1.SpectralBounds(&lower, &upper);
    if (lower <= 0.0) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
  }

  // verify monotonicity property
  AmanziGeometry::Point T1a(d);
  T1a = T1 * areas;
  for (int i = 0; i < d; i++) {
    if (T1a[i] <= 0.0) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
  }

  // add stability term D_ik T1_kl D_il
  double volume = mesh_->cell_volume(c);
  for (int i = 0; i < nrows; i++) {
    i1 = i / 2;
    k = map[i];
    area1 = mesh_->face_area(faces[k]);

    for (int j = i; j < nrows; j++) {
      i2 = j / 2;
      l = map[j];
      area2 = mesh_->face_area(faces[l]);
      M(l, k) = M(k, l) += T1(i1, i2) * area1 * area2 / volume;  // Fix (lipnikov@lanl.gov)
    }
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Simple stability term for nonsymmetric tensors.
****************************************************************** */
void MFD3D_Diffusion::StabilityScalarNonSymmetric_(int c, DenseMatrix& N, DenseMatrix& M)
{
  GrammSchmidt(N);
  CalculateStabilityScalar(M);

  int nrows = M.NumRows();
  int ncols = N.NumCols();

  // add projector ss * (I - N^T N) to matrix M
  for (int i = 0; i < nrows; i++) {  
    M(i, i) += scalar_stability_;

    for (int j = i; j < nrows; j++) {
      double s = 0.0;
      for (int k = 0; k < ncols; k++)  s += N(i, k) * N(j, k);

      s *= scalar_stability_;
      M(i, j) -= s;
      if (i - j) M(j, i) -= s;
    }
  }
}


}  // namespace WhetStone
}  // namespace Amanzi



