/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The mimetic finite difference method for diffusion.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for inner product in space of fluxes.
* Only upper triangular part of Mc = R (R^T N)^{-1} R^T is calculated.
* Here R^T N = |c| K.
* Fluxes include face areas!
****************************************************************** */
int
MFD3D_Diffusion::L2consistencyScaledArea(int c,
                                         const Tensor& K,
                                         DenseMatrix& N,
                                         DenseMatrix& Mc,
                                         bool symmetry)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_);
  Mc.Reshape(nfaces, nfaces);

  double volume = mesh_->getCellVolume(c);

  AmanziGeometry::Point v1(d_), v2(d_);
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  Tensor Kinv(K);
  Kinv.Inverse();
  Kinv.Transpose();

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    v2 = Kinv * (fm - cm);

    int i0 = (symmetry ? i : 0);
    for (int j = i0; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& fm2 = mesh_->getFaceCentroid(f);
      v1 = fm2 - cm;
      Mc(i, j) = (v1 * v2) / volume;
    }
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    for (int k = 0; k < d_; k++) N(i, k) = normal[k] * dirs[i];
  }

  return 0;
}


/* ******************************************************************
* Consistency condition for inverse of the mass matrix in the space
* of Darcy fluxes. Only the upper triangular part of matrix
* Wc = N (N^T R)^{-1} N^T is calculated. Here N^T R = |c| K.
* Flux is scaled by face area!
****************************************************************** */
int
MFD3D_Diffusion::L2consistencyInverseScaledArea(int c,
                                                const Tensor& K,
                                                DenseMatrix& R,
                                                DenseMatrix& Wc,
                                                bool symmetry)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  R.Reshape(nfaces, d_);
  Wc.Reshape(nfaces, nfaces);

  AmanziGeometry::Point v1(d_);
  double volume = mesh_->getCellVolume(c);

  Tensor Kt(K);
  Kt.Transpose();

  // Since N is scaled by K, N = N0 * K, we us tensor K in the
  // inverse L2 consistency term.
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    v1 = Kt * normal;

    int i0 = (symmetry ? i : 0);
    for (int j = i0; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->getFaceNormal(f);
      Wc(i, j) = (v1 * v2) / (dirs[i] * dirs[j] * volume);
    }
  }

  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    for (int k = 0; k < d_; k++) R(i, k) = fm[k] - cm[k];
  }

  return 0;
}


/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction.
* Only the upper triangular part of Ac is calculated.
* The degrees of freedom are at nodes.
****************************************************************** */
int
MFD3D_Diffusion::H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  auto nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  N.Reshape(nnodes, d_ + 1);
  Ac.Reshape(nnodes, nnodes);

  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);

  double volume = mesh_->getCellVolume(c);
  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  // to calculate matrix R, we use temporary matrix N
  N.PutScalar(0.0);

  int num_faces = faces.size();
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);
    double area = mesh_->getFaceArea(f);

    auto face_nodes = mesh_->getFaceNodes(f);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d_ == 2) {
        u = 0.5 * dirs[i];
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        p = mesh_->getNodeCoordinate(v);
        pnext = mesh_->getNodeCoordinate(vnext);
        pprev = mesh_->getNodeCoordinate(vprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1 ^ v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
      for (int k = 0; k < d_; k++) N(pos, k) += normal[k] * u;
    }
  }

  // calculate upper part of R K R^T / volume
  for (int i = 0; i < nnodes; i++) {
    for (int k = 0; k < d_; k++) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < nnodes; j++) {
      for (int k = 0; k < d_; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh_->getNodeCoordinate(v);
    for (int k = 0; k < d_; k++) N(i, k) = p[k] - cm[k];
    N(i, d_) = 1.0; // additional column is added to the consistency condition
  }

  return 0;
}


/* ******************************************************************
* Mass matrix in space of fluxes.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixScaledArea(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistencyScaledArea(c, Kinv, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}


/* ******************************************************************
* Mass matrix in space of fluxes for non-symmetric tensor
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixNonSymmetric(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistencyScaledArea(c, Kinv, N, M, false);
  if (ok) return ok;

  StabilityScalarNonSymmetric_(N, M);
  return 0;
}


/* ******************************************************************
* Inverse mass matrix for non-symmetric PD tensor.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseNonSymmetric(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverseScaledArea(c, K, R, W, false);
  if (ok) return ok;

  StabilityScalarNonSymmetric_(R, W);
  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm.
****************************************************************** */
int
MFD3D_Diffusion::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Stiffness matrix: the optimized algorithm.
****************************************************************** */
int
MFD3D_Diffusion::StiffnessMatrixOptimized(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  ok = StabilityOptimized_(K, N, A);
  return ok;
}


/* ******************************************************************
* Stiffness matrix: the M-matrix approach
****************************************************************** */
int
MFD3D_Diffusion::StiffnessMatrixMMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  // scaling of matrix A for numerical stability
  double s = A.Trace() / A.NumRows();
  A /= s;

  int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_TRACE;
  StabilityMMatrix_(c, N, A, objective);
  if (ok) return ok;

  A *= s;
  simplex_functional_ *= s;
  return 0;
}


/* *****************************************************************
* Low-order L2 projector.
* NOTE: we abuse the interface and return a linear polynomial.
***************************************************************** */
void
MFD3D_Diffusion::L2Cell(int c,
                        const std::vector<Polynomial>& ve,
                        const std::vector<Polynomial>& vf,
                        const Polynomial* moments,
                        Polynomial& vc)
{
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  vc.Reshape(d_, 1, true);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->getFaceCentroid(f);

    for (int k = 0; k < d_; k++) {
      double Rik = fm[k] - cm[k];
      vc(k + 1) += Rik * vf[i](0) * dirs[i];
    }
  }

  vc *= -1.0 / mesh_->getCellVolume(c);
}


/* ******************************************************************
* Divergence matrix.
****************************************************************** */
int
MFD3D_Diffusion::DivergenceMatrix(int c, DenseMatrix& A)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  A.Reshape(1, nfaces);

  for (int n = 0; n < nfaces; ++n) { A(0, n) = mesh_->getFaceArea(faces[n]) * dirs[n]; }
  return 0;
}


/* *****************************************************************
*  OTHER ROUTINES
***************************************************************** */

/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of
* fluxes. Only the upper triangular part of Wc is calculated.
****************************************************************** */
int
MFD3D_Diffusion::L2consistencyInverseDivKScaled(int c,
                                                const Tensor& K,
                                                double kmean,
                                                const AmanziGeometry::Point& kgrad,
                                                DenseMatrix& R,
                                                DenseMatrix& Wc)
{

  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  R.Reshape(nfaces, d_);
  Wc.Reshape(nfaces, nfaces);

  // calculate areas of possibly curved faces
  AmanziMesh::Double_List areas(nfaces, 0.0);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->getFaceNormal(f));
  }

  // populate matrix W_0
  AmanziGeometry::Point v1(d_), v2(d_);
  double volume = mesh_->getCellVolume(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    v1 = (K * normal) / kmean;

    for (int j = i; j < nfaces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v3 = mesh_->getFaceNormal(f);
      Wc(i, j) = (v1 * v3) / (dirs[i] * dirs[j] * volume * areas[i] * areas[j]);
    }
  }

  // populate matrix R
  const AmanziGeometry::Point& cm = mesh_->getCellCentroid(c);

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    if (d_ == 2) {
      auto nodes = mesh_->getFaceNodes(f);

      v1 = mesh_->getNodeCoordinate(nodes[0]);
      v2 = mesh_->getNodeCoordinate(nodes[1]);

      v1 -= cm;
      v2 -= cm;

      double k1 = kmean + kgrad * v1;
      double k2 = kmean + kgrad * v2;
      double km = k1 + k2;
      double tmp = areas[i] / 6;
      for (int k = 0; k < d_; k++) {
        double vm = v1[k] + v2[k];
        R(i, k) = (k1 * v1[k] + k2 * v2[k] + km * vm) * tmp;
      }
    } else {
      AMANZI_ASSERT(false);
    }
  }

  return 0;
}


/* ******************************************************************
* Inverse mass matrix in flux space via optimization, experimental.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return ok;

  StabilityScalar_(R, W);
  RescaleMassMatrixInverse_(c, W);

  return ok;
}


/* ******************************************************************
* Mass matrix for a hexahedral element, a brick for now.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseMMatrixHex(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverseScaledArea(c, K, R, W, true);
  if (ok) return ok;

  ok = StabilityMMatrixHex_(c, K, W);
  if (ok) return ok;

  return 0;
}


/* ******************************************************************
* Mass matrix for a polyhedral element via simplex method.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseMMatrix(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  // use boolean flag to populate the whole matrix
  int ok = L2consistencyInverseScaledArea(c, K, R, W, false);
  if (ok) return ok;

  // scaling of matrix W for numerical stability
  double s = W.Trace() / W.NumRows();
  W /= s;

  ok = StabilityMMatrix_(c, R, W);
  if (ok) return ok;

  W *= s;
  simplex_functional_ *= s;
  return 0;
}


/* ******************************************************************
* Inverse mass matrix via optimization, experimental.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseOptimized(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return ok;

  ok = StabilityOptimized_(K, R, W);
  RescaleMassMatrixInverse_(c, W);

  return ok;
}


/* ******************************************************************
* Inverse mass matrix via optimization, experimental.
****************************************************************** */
int
MFD3D_Diffusion::MassMatrixInverseDivKScaled(int c,
                                             const Tensor& K,
                                             double kmean,
                                             const AmanziGeometry::Point& kgrad,
                                             DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverseDivKScaled(c, K, kmean, kgrad, R, W);
  if (ok) return ok;

  StabilityScalar_(R, W);
  RescaleMassMatrixInverse_(c, W);

  return ok;
}


/* ******************************************************************
* Rescale matrix to area-weighted fluxes.
****************************************************************** */
void
MFD3D_Diffusion::RescaleMassMatrixInverse_(int c, DenseMatrix& W)
{
  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  // calculate areas of possibly curved faces
  AmanziMesh::Double_List areas(nfaces, 0.0);
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    areas[i] = norm(mesh_->getFaceNormal(f));
  }

  // back to area-weighted fluxes
  for (int i = 0; i < nfaces; i++) {
    for (int j = 0; j < nfaces; j++) W(i, j) *= areas[i] * areas[j];
  }
}


/* ******************************************************************
* A simple monotone stability term for a 2D or 3D brick element.
****************************************************************** */
int
MFD3D_Diffusion::StabilityMMatrixHex_(int c, const Tensor& K, DenseMatrix& M)
{
  int nrows = 2 * d_;

  // symmetrize the consistency matrix
  for (int i = 0; i < nrows; i++) {
    for (int j = i; j < nrows; j++) M(j, i) = M(i, j);
  }

  // create groups of quasi-parallel faces
  int map[nrows];
  for (int i = 0; i < nrows; i++) map[i] = i;

  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);

  int i1, i2, k, l;
  double s1, s2, area1, area2;
  for (int i = 0; i < d_ - 1; i++) {
    i1 = 2 * i;
    i2 = i1 + 1;

    int f = faces[i1];
    const AmanziGeometry::Point& normal1 = mesh_->getFaceNormal(f);
    area1 = mesh_->getFaceArea(f);

    s1 = 1.0;
    for (int j = i2; j < nrows; j++) {
      f = faces[j];
      const AmanziGeometry::Point& normal2 = mesh_->getFaceNormal(f);
      area2 = mesh_->getFaceArea(f);

      s2 = (normal1 * normal2) * (dirs[i] * dirs[j]) / (area1 * area2);
      if (s2 < s1) { // swap map values in positions i2 and j
        k = map[i2];
        map[i2] = j;
        map[j] = k;
        s1 = s2;
      }
    }
    // if (s1 >= 0.0) return -1;  // hex is too disturted
  }

  // define transformed tensor
  Tensor T1(d_, 2);
  AmanziGeometry::Point areas(d_);
  for (int i = 0; i < d_; i++) {
    k = map[2 * i];
    int f = faces[k];
    const AmanziGeometry::Point& normal1 = mesh_->getFaceNormal(f);
    area1 = mesh_->getFaceArea(f);
    areas[i] = area1;

    for (int j = i; j < d_; j++) {
      l = map[2 * j];
      f = faces[l];
      const AmanziGeometry::Point& normal2 = mesh_->getFaceNormal(f);
      area2 = mesh_->getFaceArea(f);

      s1 = (K * normal1) * normal2 * (dirs[k] * dirs[l]) / (area1 * area2);
      if (i - j) {
        T1(i, j) = T1(j, i) = -fabs(s1);
      } else {
        T1(i, i) = s1;
      }
    }
  }

  // verify SPD property
  if (d_ == 3) {
    double lower, upper;
    T1.SpectralBounds(&lower, &upper);
    if (lower <= 0.0) return 1;
  }

  // verify monotonicity property
  AmanziGeometry::Point T1a(d_);
  T1a = T1 * areas;
  for (int i = 0; i < d_; i++) {
    if (T1a[i] <= 0.0) return 1;
  }

  // add stability term D_ik T1_kl D_il
  double volume = mesh_->getCellVolume(c);
  for (int i = 0; i < nrows; i++) {
    i1 = i / 2;
    k = map[i];
    area1 = mesh_->getFaceArea(faces[k]);

    for (int j = i; j < nrows; j++) {
      i2 = j / 2;
      l = map[j];
      area2 = mesh_->getFaceArea(faces[l]);
      M(l, k) = M(k, l) += T1(i1, i2) * area1 * area2 / volume; // Fix (lipnikov@lanl.gov)
    }
  }
  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
