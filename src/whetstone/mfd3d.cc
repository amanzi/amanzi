/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: ara-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#include <cmath>
#include <vector>

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructors
****************************************************************** */
MFD3D::MFD3D(Teuchos::RCP<AmanziMesh::Mesh> mesh) 
{ 
  mesh_ = mesh; 
  stability_method_ = WHETSTONE_STABILITY_GENERIC;
}


/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D::DarcyMass(int cell, const Tensor& permeability,
                     Teuchos::SerialDenseMatrix<int, double>& M)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);
  int nfaces = faces.size();

  Teuchos::SerialDenseMatrix<int, double> N(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Mc(nfaces, nfaces);

  Tensor permeability_inv(permeability);
  permeability_inv.inverse();

  int ok = L2consistency(cell, permeability_inv, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, N, Mc, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D::DarcyMassInverse(int cell, const Tensor& permeability,
                            Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);
  int nfaces = faces.size();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int MFD3D::DarcyMassInverseHex(int cell, const Tensor& permeability,
                               Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);
  int nfaces = faces.size();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  int flag = StabilityMonotoneHex(cell, permeability, Wc, W);
  if (flag) StabilityScalar(cell, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* This is a debug version of the above routine for a scalar tensor
* and an orthogonal brick element.
****************************************************************** */
int MFD3D::DarcyMassInverseDiagonal(int cell, const Tensor& permeability,
                                    Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
  int nfaces = faces.size();

  W.putScalar(0.0);
  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    double area = mesh_->face_area(f);
    W(n, n) = nfaces * permeability(0, 0) * area * area / (d * volume);
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Second-generation MFD method as inlemented in RC1.
****************************************************************** */
int MFD3D::DarcyMassInverseSO(int cell, const Tensor& permeability,
                              Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);

  AmanziMesh::Entity_ID_List nodes, corner_faces;
  mesh_->cell_get_nodes(cell, &nodes);
  int nnodes = nodes.size();

  Tensor K(permeability);
  K.inverse();

  // collect all corner matrices
  std::vector<Tensor> Mv;
  std::vector<double> cwgt;

  Tensor N(d, 2), NK(d, 2), Mv_tmp(d, 2);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_cell_faces(v, cell, AmanziMesh::USED, &corner_faces);
    int nfaces = corner_faces.size();
    if (nfaces < d) {
      Errors::Message msg;
      msg << "WhetStone MFD3D: number of faces forming a corner is small.";
      Exceptions::amanzi_throw(msg);
    }

    for (int i = 0; i < d; i++) {
      int f = corner_faces[i];
      N.add_column(i, mesh_->face_normal(f));
    }
    double cwgt_tmp = fabs(N.determinant());

    N.inverse();
    NK = N * K;

    N.transpose();
    Mv_tmp = NK * N;
    Mv.push_back(Mv_tmp);

    for (int i = 0; i < d; i++) {
      int f = corner_faces[i];
      cwgt_tmp /= mesh_->face_area(f);
    }
    cwgt.push_back(cwgt_tmp);
  }

  // rescale corner weights
  double factor = 0.0;
  for (int n = 0; n < nnodes; n++) factor += cwgt[n];
  factor = mesh_->cell_volume(cell) / factor;

  for (int n = 0; n < nnodes; n++) cwgt[n] *= factor;

  // assemble corner matrices
  W.putScalar(0.0);
  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    mesh_->node_get_cell_faces(v, cell, AmanziMesh::USED, &corner_faces);

    Tensor& Mv_tmp = Mv[n];
    for (int i = 0; i < d; i++) {
      int k = FindPosition(corner_faces[i], faces);
      for (int j = i; j < d; j++) {
        int l = FindPosition(corner_faces[j], faces);
        W(k, l) += Mv_tmp(i, j) * cwgt[n] * fdirs[k] * fdirs[l];
        W(l, k) = W(k, l);
      }
    }
  }
 
  // invert matrix W
  Teuchos::LAPACK<int, double> lapack;
  int info, size = W.numRows();

  int ipiv[size];
  double work[size];

  lapack.GETRF(size, size, W.values(), size, ipiv, &info);
  lapack.GETRI(size, W.values(), size, ipiv, work, size, &info);
  if (info != 0) {
    Errors::Message msg;
    msg << "WhetStone MFD3D: support operator generated bad elemental mass matrix.";
    Exceptions::amanzi_throw(msg);
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int MFD3D::DarcyMassInverseOptimized(int cell, const Tensor& permeability,
                                     Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);
  int nfaces = faces.size();

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2consistencyInverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  ok = StabilityOptimized(permeability, R, Wc, W);
  return ok;
}


/* ******************************************************************
* Consistency condition for inner product in space of Darcy fluxes. 
* Only upper triangular part of Mc is calculated.
* Darcy flux is scaled by area!
****************************************************************** */
int MFD3D::L2consistency(int cell, const Tensor& T,
                         Teuchos::SerialDenseMatrix<int, double>& N,
                         Teuchos::SerialDenseMatrix<int, double>& Mc)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  mesh_->cell_get_faces_and_dirs(cell, &faces, &fdirs);

  int num_faces = faces.size();
  if (num_faces != N.numRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);

  AmanziGeometry::Point v1(d), v2(d);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  Tensor Tinv(T);
  Tinv.inverse();

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    v2 = Tinv * (fm - cm);

    for (int j = i; j < num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
      v1 = fm - cm;
      Mc(i, j) = (v1 * v2) / volume;
    }
  }

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    for (int k = 0; k < d; k++) N(i, k) = normal[k];
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of Darcy 
* fluxes. Only the upper triangular part of Wc is calculated.
* Darcy flux is scaled by area!
****************************************************************** */
int MFD3D::L2consistencyInverse(int cell, const Tensor& T,
                                Teuchos::SerialDenseMatrix<int, double>& R,
                                Teuchos::SerialDenseMatrix<int, double>& Wc)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.numRows()) return num_faces;  // matrix was not reshaped

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
int MFD3D::H1consistency(int cell, const Tensor& T,
                         Teuchos::SerialDenseMatrix<int, double>& N,
                         Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(cell, &nodes);
  int num_nodes = nodes.size();
  if (num_nodes != N.numRows()) return num_nodes;  // matrix was not reshaped

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  /* to calculate matrix R, we use temporary matrix N */
  N = 0;

  int num_faces = faces.size();
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    AmanziMesh::Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int jnext = (j + 1) % num_face_nodes;
      int jprev = (j + num_face_nodes - 1) % num_face_nodes;

      int v = face_nodes[j];
      int vnext = face_nodes[jnext];
      int vprev = face_nodes[jprev];

      mesh_->node_get_coordinates(v, &p);
      mesh_->node_get_coordinates(vnext, &pnext);
      mesh_->node_get_coordinates(vprev, &pprev);

      v1 = pprev - pnext;
      v2 = p - fm;
      v3 = v1^v2;
      double u = dirs[i] * norm(v3) / (4 * area);

      int pos = FindPosition(v, nodes);
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
* Consistency condition for stifness matrix in geomechanics. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D::H1consistencyElasticity(int cell, const Tensor& T,
                                   Teuchos::SerialDenseMatrix<int, double>& N,
                                   Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  return WHETSTONE_ELEMENTAL_MATRIX_OK;  // (lipnikov@lanl.gov)
}


/* ******************************************************************
* Calculate stability factor using matrix and optional scaling.
****************************************************************** */
double MFD3D::CalculateStabilityScalar(Teuchos::SerialDenseMatrix<int, double>& Mc)
{
  int nrows = Mc.numRows();

  scalar_stability_ = 0.0;
  for (int i = 0; i < nrows; i++) scalar_stability_ += Mc(i, i);
  scalar_stability_ /= nrows;

  if (stability_method_ == WHETSTONE_STABILITY_GENERIC_SCALED) {
    scalar_stability_ *= scaling_factor_;
  }

  return scalar_stability_;
}


/* ******************************************************************
* Set up positive scaling factor for a scalar stability term.
* Warning: Ignores silently negative factors.
****************************************************************** */
void MFD3D::ModifyStabilityScalingFactor(double factor)
{
  if (factor > 0.0) {
    stability_method_ = WHETSTONE_STABILITY_GENERIC_SCALED;
    scaling_factor_ = factor;
  }
}


/* ******************************************************************
* Simplest stability term is added to the consistency term. 
****************************************************************** */
void MFD3D::StabilityScalar(int cell,
                            Teuchos::SerialDenseMatrix<int, double>& N,
                            Teuchos::SerialDenseMatrix<int, double>& Mc,
                            Teuchos::SerialDenseMatrix<int, double>& M)
{
  GrammSchmidt(N);
  CalculateStabilityScalar(Mc);

  int nrows = Mc.numRows();
  int ncols = N.numCols();

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(i, j) = Mc(i, j);
  }

  for (int i = 0; i < nrows; i++) {  // add projector ss * (I - N^T N) to matrix M
    M(i, i) += scalar_stability_;

    for (int j = i; j < nrows; j++) {
      double s = 0.0;
      for (int k = 0; k < ncols; k++)  s += N(i, k) * N(j, k);
      M(i, j) -= s * scalar_stability_;
    }
  }

  for (int i = 0; i < nrows; i++) {  // symmetrization
    for (int j = i+1; j < nrows; j++) M(j, i) = M(i, j);
  }
}


/* ******************************************************************
* A simple monotone stability term for a 2D or 3D brick element. 
****************************************************************** */
int MFD3D::StabilityMonotoneHex(int cell, const Tensor& T,
                                Teuchos::SerialDenseMatrix<int, double>& Mc,
                                Teuchos::SerialDenseMatrix<int, double>& M)
{
  int d = mesh_->space_dimension();
  int nrows = 2 * d;

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(j, i) = M(i, j) = Mc(i, j);
  }

  // create groups of quasi-parallel faces
  int map[nrows];
  for (int i = 0; i < nrows; i++) map[i] = i;

  AmanziMesh::Entity_ID_List faces;
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
  for (int i = 0; i < d; i++) {
    k = map[2*i];
    int f = faces[k];
    const AmanziGeometry::Point& normal1 = mesh_->face_normal(f);
    area1 = mesh_->face_area(f);

    for (int j = i; j < d; j++) {
      l = map[2*j];
      f = faces[l];
      const AmanziGeometry::Point& normal2 = mesh_->face_normal(f);
      area2 = mesh_->face_area(f);

      s1 = (T * normal1) * normal2 * (dirs[k] * dirs[l]) / (area1 * area2);
      if (i-j)
        T1(i, j) = T1(j, i) = -fabs(s1);
      else
        T1(i, i) = s1;
    }
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
  return 0;
}


/* ******************************************************************
* A simple optimization procedure that must return a diagonal mass
* matrix for a 2D and 3D orthogonal cells and diagonal tensors. 
* The algorithm minimizes off-diagonal entries in the mass matrix.
* WARNING: the routine is used for inverse of mass matrix only.
****************************************************************** */
int MFD3D::StabilityOptimized(const Tensor& T,
                              Teuchos::SerialDenseMatrix<int, double>& N,
                              Teuchos::SerialDenseMatrix<int, double>& Mc,
                              Teuchos::SerialDenseMatrix<int, double>& M)
{
  int d = mesh_->space_dimension();
  int nrows = N.numRows();
  int ncols = N.numCols();

  // find correct scaling of a stability term
  double lower, upper, eigmin = Mc(0, 0);
  // T.spectral_bounds(&lower, &upper);
  for (int k = 1; k < nrows; k++) eigmin = std::min<double>(eigmin, Mc(k, k));

  // find null space of N^T
  Teuchos::SerialDenseMatrix<int, double> U(nrows, nrows);
  int info, size = 5 * d + 3 * nrows;
  double V, S[nrows], work[size];

  Teuchos::LAPACK<int, double> lapack;
  lapack.GESVD('A', 'N', nrows, ncols, N.values(), nrows,  // N = u s v
               S, U.values(), nrows, &V, 1, work, size, 
               NULL, &info);
  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // calculate vectors C and C0
  int mrows = nrows * (nrows - 1) / 2;
  int mcols = nrows - ncols;
  int nparam = (mcols + 1) * mcols / 2;
  Teuchos::SerialDenseMatrix<int, double> C(mrows, nparam);
  Teuchos::SerialDenseVector<int, double> F(mrows);

  int m, n = 0;
  for (int k = ncols; k < nrows; k++) {
    m = 0;  // calculate off-diagonal entries of M_kk = U_k * U_k^T
    for (int i = 0; i < nrows; i++) 
      for (int j = i+1; j < nrows; j++) C(m++, n) = U(i, k) * U(j, k);
    n++; 
  }

  for (int k = ncols; k < nrows; k++) {
    for (int l = k+1; l < nrows; l++) {
      m = 0;  // calculate off-diagonal entries of M_kk + M_ll - M_kl - M_lk 
      for (int i = 0; i < nrows; i++) { 
        for (int j = i+1; j < nrows; j++) {
          C(m, n) = C(m, k-ncols) + C(m, l-ncols) - U(i, k) * U(j, l) - U(i, l) * U(j, k);
          m++;
        }
      }
      n++;
    }
  }

  m = 0;
  for (int i = 0; i < nrows; i++) { 
    for (int j = i+1; j < nrows; j++) F(m++) = -Mc(i, j);
  }

  // Form a linear system for parameters
  Teuchos::SerialDenseMatrix<int, double> A(nparam, nparam);
  Teuchos::SerialDenseVector<int, double> G(nparam);

  A.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, C, C, 0.0);
  G.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, C, F, 0.0);

  // Find parameters
  lapack.POSV('U', nparam, 1, A.values(), nparam, G.values(), nparam, &info);
  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // project solution on the positive quadrant and convert to matrix
  Teuchos::SerialDenseMatrix<int, double> P(mcols, mcols);

  int status = WHETSTONE_ELEMENTAL_MATRIX_OK;
  for (int loop = 0; loop < 3; loop++) {
    if (loop == 1) {   
      for (int i = 0; i < mcols; i++) G(i) = std::max<double>(G(i), 0.0);
      status = WHETSTONE_ELEMENTAL_MATRIX_PASSED;
    } else if (loop == 2) {
      for (int i = mcols; i < nparam; i++) G(i) = std::max<double>(G(i), 0.0);
      status = WHETSTONE_ELEMENTAL_MATRIX_PASSED;
    }

    for (int k = 0; k < mcols; k++) P(k, k) = G(k);

    n = mcols;
    for (int k = 0; k < mcols; k++) {
      for (int l = k+1; l < mcols; l++) {
        P(k, k) += G(n);
        P(l, l) += G(n);
        P(l, k) = P(k, l) = -G(n);
        n++;  
      }
    }

    // check SPD property (we use allocated memory)
    Teuchos::SerialDenseMatrix<int, double> Ptmp(P);
    lapack.SYEV('N', 'U', mcols, Ptmp.values(), mcols, S, work, size, &info); 
    if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

    if (S[0] > eigmin) {
      break;
    } else if (loop == 2) {
      for (int k = 0; k < mcols; k++) if (P(k, k) == 0.0) P(k, k) = eigmin;
    }
  }

  // add stability term U G U^T
  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(i, j) = Mc(i, j);
  }

  Teuchos::SerialDenseMatrix<int, double> UP(nrows, mcols);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < mcols; j++) {
      double& entry = UP(i, j);
      for (int k = 0; k < mcols; k++) entry += U(i, k+ncols) * P(k, j);
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) {
      double& entry = M(i, j);
      for (int k = 0; k < mcols; k++) entry += UP(i, k) * U(j, k+ncols);
      M(j, i) = M(i, j); 
    }
  }

  return status;
}


/* ******************************************************************
* Conventional Gramm-Schimdt orthogonalization of colums of matrix N. 
****************************************************************** */
void MFD3D::GrammSchmidt(Teuchos::SerialDenseMatrix<int, double>& N)
{
  int nrows = N.numRows();
  int ncols = N.numCols();

  int i, j, k;
  for (i = 0; i < ncols; i++) {
    double l22 = 0.0;
    for (k = 0; k < nrows; k++) l22 += N(k, i) * N(k, i);

    l22 = 1.0 / sqrt(l22);
    for (k = 0; k < nrows; k++) N(k, i) *= l22;

    for (j = i+1; j < ncols; j++) {
      double s = 0.0;
      for (k = 0; k < nrows; k++) s += N(k, i) * N(k, j);
      for (k = 0; k < nrows; k++) N(k, j) -= s * N(k, i);  // orthogonolize i and j
    }
  }
}


/* ******************************************************************
* Extension of Mesh API. 
****************************************************************** */
int MFD3D::cell_get_face_adj_cell(const int cell, const int face)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(face, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  if (ncells == 2) {
    int c2 = cells[0];
    if (cell == c2) c2 = cells[1];
    return c2;
  }
  return -1;
}


/* ******************************************************************
* Returns position of the number v in the list of nodes.  
****************************************************************** */
int MFD3D::FindPosition(int v, AmanziMesh::Entity_ID_List nodes)
{
  for (int i = 0; i < nodes.size(); i++) {
    if (nodes[i] == v) return i;
  }
  return -1;
}

}  // namespace WhetStone
}  // namespace Amanzi



