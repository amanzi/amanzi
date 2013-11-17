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
#include "Point.hh"
#include "errors.hh"

#include "mfd3d_diffusion.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for inner product in space of fluxes. 
* Only upper triangular part of Mc is calculated.
* Darcy flux is scaled by the area!
****************************************************************** */
int MFD3D_Diffusion::L2consistency(int cell, const Tensor& T,
                                   DenseMatrix& N, DenseMatrix& Mc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int nfaces = faces.size();
  if (nfaces != N.NumRows()) return nfaces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);

  AmanziGeometry::Point v1(d), v2(d);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    v2 = Tinv * (fm - cm);

    for (int j = i; j < nfaces; j++) {
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
* Consistency condition for inverse of mass matrix in space of
* fluxes. Only the upper triangular part of Wc is calculated.
* Darcy flux is scaled by the area!
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverse(int cell, const Tensor& T,
                                          DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.NumRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(cell);

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


/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_Diffusion::H1consistency(int cell, const Tensor& T,
                                   DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(cell, &nodes);
  int num_nodes = nodes.size();
  if (num_nodes != N.NumRows()) return num_nodes;  // matrix was not reshaped

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);
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

  for (int i = 0; i < num_nodes; i++) {  // calculate R T R^T / volume
    for (int k = 0; k < d; k++) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j = i; j < num_nodes; j++) {
      for (int k = 0; k < d; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);
  for (int i = 0; i < num_nodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    for (int k = 0; k < d; k++) N(i, k) = p[k] - cm[k];
    N(i, d) = 1;  // additional colum is added to the consistency condition
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix: the standard algorithm
****************************************************************** */
int MFD3D_Diffusion::MassMatrix(int cell, const Tensor& permeability, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nfaces = M.NumRows();

  DenseMatrix N(nfaces, d);
  DenseMatrix Mc(nfaces, nfaces);

  Tensor permeability_inv(permeability);
  permeability_inv.Inverse();

  int ok = L2consistency(cell, permeability_inv, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, N, Mc, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy inverse mass matrix: the standard algorithm
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverse(int cell, const Tensor& permeability, 
                                       DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);
  DenseMatrix Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy stiffness matrix: the standard algorithm.
****************************************************************** */
int MFD3D_Diffusion::StiffnessMatrix(int cell, const Tensor& permeability, 
                                     DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nnodes = A.NumRows();

  DenseMatrix N(nnodes, d + 1);
  DenseMatrix Ac(nnodes, nnodes);

  int ok = H1consistency(cell, permeability, N, Ac);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, N, Ac, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy stiffness matrix: the M-matrix approach
****************************************************************** */
int MFD3D_Diffusion::StiffnessMatrixMMatrix(int cell, const Tensor& permeability, 
                                            DenseMatrix& A)
{
  int d = mesh_->space_dimension();
  int nnodes = A.NumRows();

  DenseMatrix N(nnodes, d + 1);
  DenseMatrix Ac(nnodes, nnodes);

  int ok = H1consistency(cell, permeability, N, Ac);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityMMatrix_(cell, N, Ac, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* *****************************************************************
* Recover gradient from solution, which is either the edge-based 
* fluxes of node-based pressures. The algorithm is common if both
* N and R are used. Here we use simplified versions.
***************************************************************** */
int MFD3D_Diffusion::RecoverGradient_MassMatrix(int cell,
                                                const std::vector<double>& solution, 
                                                AmanziGeometry::Point& gradient)
{
  int d = mesh_->space_dimension();
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
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

  gradient *= -1.0 / mesh_->cell_volume(cell);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* *****************************************************************
* Recover gradient from solution, which is either the edge-based 
* fluxes of node-based pressures. The algorithm is common if both
* N and R are used. Here we use simplified versions.
***************************************************************** */
int MFD3D_Diffusion::RecoverGradient_StiffnessMatrix(int cell,
                                                     const std::vector<double>& solution, 
                                                     AmanziGeometry::Point& gradient)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(cell, &nodes);
  int num_nodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
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
  gradient *= 1.0 / mesh_->cell_volume(cell);

  return 0;
}


/* *****************************************************************
*  OTHER ROUTINES
***************************************************************** */

/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of Darcy 
* fluxes. Only the upper triangular part of Wc is calculated.
****************************************************************** */
int MFD3D_Diffusion::L2consistencyInverseScaled(int cell, const Tensor& T,
                                                DenseMatrix& R, DenseMatrix& Wc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

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
  double volume = mesh_->cell_volume(cell);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    v1 = T * normal;

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

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
* Darcy inverse mass matrix via optimization, experimental.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseScaled(int cell, const Tensor& permeability,
                                             DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);
  DenseMatrix Wc(nfaces, nfaces);

  int ok = L2consistencyInverseScaled(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
 
  StabilityScalar(cell, R, Wc, W);
  RescaleMassMatrixInverse_(cell, W);

  return ok;
}


/* ******************************************************************
* Darcy mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseMMatrixHex(
    int cell, const Tensor& permeability, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);
  DenseMatrix Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityMMatrixHex_(cell, permeability, Wc, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
}


/* ******************************************************************
* Darcy mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseMMatrix(
    int cell, const Tensor& permeability, DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);
  DenseMatrix Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityMMatrix_(cell, R, Wc, W);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy inverse mass matrix via optimziation.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseOptimized(int cell, const Tensor& permeability,
                                                DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);
  DenseMatrix Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(permeability, R, Wc, W);
  return ok;
}


/* ******************************************************************
* Darcy inverse mass matrix via optimization, experimental.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseOptimizedScaled(int cell, const Tensor& permeability,
                                                      DenseMatrix& W)
{
  int d = mesh_->space_dimension();
  int nfaces = W.NumRows();

  DenseMatrix R(nfaces, d);
  DenseMatrix Wc(nfaces, nfaces);

  int ok = L2consistencyInverseScaled(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;
 
  ok = StabilityOptimized(permeability, R, Wc, W);
  RescaleMassMatrixInverse_(cell, W);

  return ok;
}


/* ******************************************************************
* Rescale matrix to area-weighted fluxes.
****************************************************************** */
void MFD3D_Diffusion::RescaleMassMatrixInverse_(int cell, DenseMatrix& W)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
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
int MFD3D_Diffusion::StabilityMMatrixHex_(int cell, const Tensor& T,
                                          DenseMatrix& Mc, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = 2 * d;

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(j, i) = M(i, j) = Mc(i, j);
  }

  // create groups of quasi-parallel faces
  int map[nrows];
  for (int i = 0; i < nrows; i++) map[i] = i;

  Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

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

      s1 = (T * normal1) * normal2 * (dirs[k] * dirs[l]) / (area1 * area2);
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
  double volume = mesh_->cell_volume(cell);
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
* A wrapper for the simplex method that finds monotone parameters. 
****************************************************************** */
int MFD3D_Diffusion::StabilityMMatrix_(
    int cell, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = N.NumRows();
  int ncols = N.NumCols();

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(j, i) = M(i, j) = Mc(i, j);
  }

  // compute null space
  int mcols = nrows - ncols;
  DenseMatrix D(nrows, mcols);
  int ierr = N.NullSpace(D);

  // populate the tableau
  int m1(0), m2(nrows), m12, n12, mx, nx, ir;
  double tmp;

  m12 = nrows * (nrows + 1) / 2;
  mx = mcols * mcols;

  // Simplex method requires one auxiliary row in the tableau.
  DenseMatrix T(m12 + 2, mx + 1);
  T.PutScalar(0.0);

  // first condition M_ij < 0
  n12 = m12 - nrows;
  for (int i = 0; i < nrows; i++) {
    for (int j = i + 1; j < nrows; j++) {
      double b = M(i, j);
      if (b < 0.0) {
        ir = ++m1;
      } else {
        ir = n12--;
        m2++;
      }
      T(ir, 0) = b;
      
      nx = 0;
      for (int k = 0; k < mcols; k++) {
        tmp = D(i, k) * D(j, k);
        T(ir, ++nx) = tmp;
        for (int l = k + 1; l < mcols; l++) {
          tmp = D(i, k) * D(j, l) + D(j, k) * D(i, l);
          T(ir, ++nx) = tmp;
          T(ir, ++nx) = -tmp;
        }
      }

      if (b < 0.0) {
        for (int k = 0; k < mx + 1; k++) T(ir, k) *= -1;
      }
    }
  }

  // second condition sum_j M_ij > 0
  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) {
      nx = 0;
      for (int k = 0; k < mcols; k++) {
        tmp = D(i, k) * D(j, k);
        nx++;
        T(m12 - i, nx) -= tmp;
        if (i != j) T(m12 - j, nx) -= tmp;

        for (int l = k + 1; l < mcols; l++) {
          tmp = D(i, k) * D(j, l) + D(j, k) * D(i, l);
          nx++;
          T(m12 - i, nx) -= tmp;
          if (i != j) T(m12 - j, nx) -= tmp;

          nx++;
          T(m12 - i, nx) += tmp;
          if (i != j) T(m12 - j, nx) += tmp;
        }
      }
    }
  }

  // objective functional
  n12 = m12 - nrows + 1;
  for (int k = 0; k < mx + 1; k++) {
    double q1 = 0.0;
    for (int i = n12; i < m12; i++) q1 += T(i, k);
    T(0, k) = -q1;
  }

  // find a feasible basic solution
  int izrow[mx + 1], iypos[m12 + 1], itrs;
  itrs = SimplexFindFeasibleSolution_(T, m1, m2, 0, izrow, iypos);
  if (itrs < 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;
  simplex_functional_ = T(0,0); 
  simplex_num_itrs_ = itrs; 

  double u[mx];
  for (int i = 0; i < mx; i++) u[i] = 0.0;
  for (int i = 1; i < m12 + 1; i++) {
    int k = iypos[i] - 1;
    if (k < mx) u[k] = T(i,0);
  }

  // add matrix D' U D
  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) { 
      nx = 0;
      for (int k = 0; k < mcols; k++) {
        M(i, j) += D(i, k) * D(j, k) * u[nx];
        nx++;
        for (int l = k + 1; l < mcols; l++) {
          double tmp = D(i, k) * D(j, l) + D(j, k) * D(i, l);
          M(i, j) += tmp * (u[nx] - u[nx + 1]);
          nx += 2;
        }
      }
      M(j, i) = M(i, j);
    }
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* A simplex method for fining monotone parameters. 
* We assume that m3 = 0; otherwise, routine MaxRowValue() has 
* to be modified by looping over columns in array l1.
****************************************************************** */
int MFD3D_Diffusion::SimplexFindFeasibleSolution_(
    DenseMatrix& T, int m1, int m2, int m3, int* izrow, int* iypos)
{
  int m = m1 + m2 + m3;     // Number of constraints.
  int n = T.NumCols() - 1;  // Number of unknowns.
  
  for (int k = 0; k < n + 1; k++) {
    double q1 = 0.0;
    for (int i = m1 + 1; i < m + 1; i++) q1 += T(i, k);
    T(m + 1, k) = -q1;
  }

  // work memory
  int ip, kp, itr_max = WHETSTONE_SIMPLEX_MAX_ITERATIONS * n;
  int itr1(0), itr2(0), nl1(n), l1[n + 1];

  for (int k = 0; k < n + 1; k++) l1[k] = izrow[k] = k;
  for (int i = 0; i < m + 1; i++) iypos[i] = n + i;

  // Start of phase I.
  if (m2 + m3 > 0) {
    int flag(0), l3[m2];
    for (int i = 0; i < m2; i++) l3[i] = 1;

    for (int itr = 0; itr < itr_max; itr++) {
      // find maximum coeffient of the auxiliary functional
      double vmax;
      T.MaxRowValue(m + 1, 1, n, &kp, &vmax); 

      // feasible solution does not exist 
      if (vmax < WHETSTONE_SIMPLEX_TOLERANCE && 
          T(m + 1, 0) < -WHETSTONE_SIMPLEX_TOLERANCE) 
          return WHETSTONE_SIMPLEX_NO_FEASIBLE_SET;

      // feasible solution has been found
      if (vmax < WHETSTONE_SIMPLEX_TOLERANCE && 
          fabs(T(m + 1, 0)) < WHETSTONE_SIMPLEX_TOLERANCE) {
        /*
        for (int ip = m1 + m2 + 1; ip < m + 1; ip++) {
          if (iypos[ip] == ip + n) {
            // Found an artificial variable for an equality constraint.
            T.MaxRowMagnitude(ip, 1, n, &kp, &vmax);
            if (vmax > WHETSTONE_SIMPLEX_TOLERANCE) goto one;
          }
        }
        */

        for (int i = m1 + 1; i < m + 1; i++) {
          if (l3[i - m1 - 1] == 1) {
            for (int k = 0; k < n + 1; k++) T(i, k) *= -1; 
          }
        }
        itr1 = itr;
        flag = 1;
        break;
      }

      // locate a pivot element in column kp (skipping degeneracy)
      SimplexPivotElement_(T, kp, &ip);
      if (ip == 0) return WHETSTONE_SIMPLEX_UNBOUNDED_PROBLEM;

      // Exchange left and right-hand variables
      SimplexExchangeVariables_(T, kp, ip);

      // Exchanged out an artificial variable foranequality constraint.
      // Make sure it stays out by removing it from the l1 list.
      if (iypos[ip] >= n + m1 + m2 + 1) {
        for (int k = 1; k <= nl1; k++) {
          if (l1[k] == kp) { 
            --nl1;
            for (int i = k; i <= nl1; i++) l1[i] = l1[i + 1];
            break;
          }
        }
      } else {
      // Exchanged out an m2 type constraint for the first time. 
      // Correct sign of the pivot column and the implicit artificial variable.
        int kh = iypos[ip] - m1 - n - 1;
        if (kh >= 0 && l3[kh] == 1) {
          l3[kh] = 0;
          T(m + 1, kp) += 1.0;
          for (int i = 0; i < m + 2; i++) T(i, kp) *= -1;
        }
      }
      // Update lists of left-hand and right-hand variables.
      int is = izrow[kp];
      izrow[kp] = iypos[ip];
      iypos[ip] = is;
    }
    if (flag == 0) return WHETSTONE_SIMPLEX_NO_CONVERGENCE;
  }

  // Start of phase II.
  for (int itr = 0; itr < itr_max; itr++) {
    double vmax;
    T.MaxRowValue(0, 1, n, &kp, &vmax);

    // solution has been found
    if (vmax < WHETSTONE_SIMPLEX_TOLERANCE) {
      itr2 = itr;
      break;
    }

    // Locate a pivot element.
    SimplexPivotElement_(T, kp, &ip);
    if (ip == 0) return WHETSTONE_SIMPLEX_UNBOUNDED_PROBLEM;

    // Exchange a left-hand and a right-hand variables.
    SimplexExchangeVariables_(T, kp, ip);

    int is = izrow[kp];
    izrow[kp] = iypos[ip];
    iypos[ip] = is;
  }

  return itr1 + itr2;
}


/* ******************************************************************
* Locates a pivot elements taking degeneracy into account.
****************************************************************** */
void MFD3D_Diffusion::SimplexPivotElement_(DenseMatrix& T, int kp, int* ip)
{
  int m = T.NumRows() - 2;
  int n = T.NumCols() - 1;
  double qmin, q, qp, q0, tmp;

  *ip = 0;
  for (int i = 1; i < m + 1; i++) {
    tmp = T(i, kp);
    if (tmp < -WHETSTONE_SIMPLEX_TOLERANCE) {
      q = -T(i, 0) / tmp;

      if (*ip == 0) {
        *ip = i;
        qmin = q;
      } else if (q < qmin) {
        *ip = i;
        qmin = q;
      } else if (q == qmin) {  // we have a degeneracy.
        double tmp0 = T(*ip, kp);
        for (int k = 1; k <= n; k++) {
          qp = -T(*ip, k) / tmp0;
          q0 = -T(i, k) / tmp;
          if (q0 != qp) break;
        }
        if (q0 < qp) *ip = i;
      }
    }
  }
}


/* ******************************************************************
* Exchanges lenf-hand and rihgt-hand variables.
****************************************************************** */
void MFD3D_Diffusion::SimplexExchangeVariables_(DenseMatrix& T, int kp, int ip)
{
  int m = T.NumRows() - 2;
  int n = T.NumCols() - 1;

  double tmp = 1.0 / T(ip, kp);
  for (int i = 0; i < m + 2; i++) {
    if (i != ip) {
      T(i, kp) *= tmp;
      for (int k = 0; k < n + 1; k++) {
        if (k != kp) T(i, k) -= T(ip, k) * T(i, kp);
      }
    }
  }
  for (int k = 0; k < n + 1; k++) {
    if (k != kp) T(ip, k) *= -tmp;
  }
  T(ip, kp) = tmp;
}

}  // namespace WhetStone
}  // namespace Amanzi



