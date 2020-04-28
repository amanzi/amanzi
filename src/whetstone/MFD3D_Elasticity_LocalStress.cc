/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method for elasticity: local stress approach
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_Elasticity.hh"
#include "Tensor.hh"
#include "WhetStoneMeshUtils.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrix_LocalStress(
   int v, const std::vector<Tensor>& T, DenseMatrix& A, DenseMatrix& B)
{
  DenseMatrix M, D, S1, S2;
  LocalStressMatrices_(v, T, M, D, S1, S2);

  DenseMatrix DT(D.NumCols(), D.NumRows());
  DT.Transpose(D);
  M.InverseMoorePenrose();

  auto Q11 = DT * M * D;
  auto Q12 = DT * M * S1;
  auto Q22 = S2 * M * S1;
  auto Q21 = S2 * M * D;
  Q22.InverseMoorePenrose();

  A = Q11 - Q12 * Q22 * Q21;
  B = (DT - Q12 * Q22 * S2) * M;
  
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix in the space of stresses.
****************************************************************** */
void MFD3D_Elasticity::LocalStressMatrices_(
   int v, const std::vector<Tensor>& T,
   DenseMatrix& M, DenseMatrix& D, DenseMatrix& S1, DenseMatrix& S2)
{
  Entity_ID_List cells, faces, cfaces, vcfaces;
  std::vector<int> cdirs;

  mesh_->node_get_faces(v, AmanziMesh::Parallel_type::ALL, &faces);
  mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
  int nfaces = faces.size();
  int ncells = cells.size();
  AMANZI_ASSERT(ncells > 0);

  int nd = d_ * (d_ + 1) / 2;  // symmetric tensors
  int nk = d_ * (d_ - 1) / 2;  // skew-symmetric tensors
  M.Reshape(d_ * nfaces, d_ * nfaces);
  D.Reshape(d_ * nfaces, d_ * ncells);
  S1.Reshape(d_ * nfaces, nk);
  S2.Reshape(nk, d_ * nfaces);

  M.PutScalar(0.0);
  D.PutScalar(0.0);
  S1.PutScalar(0.0);
  S2.PutScalar(0.0);

  // basis of stresses (symmetric stresses + non-symmetric)
  std::vector<Tensor> vE;

  for (int i = 0; i < d_; ++i) {
    Tensor E(d_, 2);
    E(i, i) = 1.0;
    vE.push_back(E);
  }

  for (int i = 0; i < d_; ++i) {
    for (int j = i + 1; j < d_; ++j) {
      Tensor E(d_, 2);
      E(i, j) = E(j, i) = 1.0;
      vE.push_back(E);
    }
  }

  for (int i = 0; i < d_; ++i) {
    for (int j = i + 1; j < d_; ++j) {
      Tensor E(d_, 2);
      E(i, j) = 1.0;
      E(j, i) =-1.0;
      vE.push_back(E);
    }
  }

  // node
  AmanziGeometry::Point xv(d_);
  mesh_->node_get_coordinates(v, &xv);

  // volume of dual cell
  double utot, volume_dual(0.0);
  std::vector<double> u(ncells, 0.0), volume_corner(ncells, 0.0);

  for (int n = 0; n < ncells; ++n) {
    int c = cells[n];
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    mesh_->node_get_cell_faces(v, c, Parallel_type::ALL, &vcfaces);
    int nvc = vcfaces.size();

    for (int i = 0; i < nvc; i++) {
      int f = vcfaces[i];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      volume_corner[n] += norm((xf - xc) ^ (xf - xv)) / 2;
    }
    volume_dual += volume_corner[n];
  }

  // main loop over cells around the given node
  for (int n = 0; n < ncells; ++n) {
    int c = cells[n];
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    mesh_->cell_get_faces_and_dirs(c, &cfaces, &cdirs);
    mesh_->node_get_cell_faces(v, c, Parallel_type::ALL, &vcfaces);
    int nvc = vcfaces.size();
    AMANZI_ASSERT(nvc == d_);

    // -- generate auxiliary corner matrix for one component
    int mx = nvc * d_, nx = d_ * d_;
    DenseMatrix Mcorner(mx, mx), Ncorner(mx, nx), Rcorner(mx, nx), Tcorner(mx, nx);
    DenseVector Dcorner(nvc);

    Tcorner.PutScalar(0.0);

    for (int i = 0; i < nvc; i++) {
      int f = vcfaces[i];
      int m = std::distance(cfaces.begin(), std::find(cfaces.begin(), cfaces.end(), f));

      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      double facet_area = area / 2;  // FIXME
      Dcorner(i) = cdirs[m] * facet_area;

      for (int k = 0; k < d_; ++k) {
        for (int l = 0; l < d_; ++l) {
          Tcorner(nvc * k + l, d_ * l + i) = normal[k] / area;
        }
      }

      for (int j = 0; j < d_ * d_; ++j) {
        auto conormal = (T[n] * vE[j]) * (normal / area);
        // auto dx = vE[j] * ((2 * xf + xv) / 3 - xc);
        auto dx = vE[j] * (xf - xc);
        // auto dx = vE[j] * ((xf + xv) / 2 - xc);

        for (int k = 0; k < d_; ++k) {
          Ncorner(nvc * k + i, j) = conormal[k];
          Rcorner(nvc * k + i, j) = cdirs[m] * dx[k] * facet_area;
        }
      }
    }

    // dual-volume scheme
    Tcorner.InverseMoorePenrose();
    DenseMatrix TcornerT(Tcorner);
    TcornerT.Transpose();
    auto Ycorner = Tcorner * TcornerT * Ncorner;

    // original scheme
    DenseMatrix Ncorner_copy(Ncorner);
    Ncorner.InverseMoorePenrose();
    Mcorner = Rcorner * Ncorner;

    // scheme based on rescaling
    DenseVector one(mx);
    one.PutScalar(0.0);
    for (int r = 0; r < mx; ++r)
      for (int j = 0; j < 1; ++j) one(r) += Ncorner_copy(r, j);

    for (int i = 0; i < mx; ++i)
      for (int j = 0; j < mx; ++j) u[n] += Mcorner(i, j) * one(i) * one(j);
    utot += u[n];

    // support-operator scheme
    // Mcorner = (Tcorner * TcornerT) * volume_corner[n];

    // assemble mass matrices
    for (int i = 0; i < nvc; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));

      for (int j = 0; j < nvc; ++j) {
        int l = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[j]));

        for (int s = 0; s < d_; ++s) { 
          for (int r = 0; r < d_; ++r) { 
            M(d_ * k + s, d_ * l + r) += Mcorner(nvc * s + i, nvc * r + j);
          }
        }
      }
    }

    // assemble divergence matrices
    for (int i = 0; i < nvc; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));
      for (int s = 0; s < d_; ++s) { 
        D(d_ * k + s, d_ * n + s) = Dcorner(i);
      }
    }

    // assemble rotation matrices
    for (int m = 0; m < nk; ++m) {
      for (int i = 0; i < nvc; i++) {
        int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));
        for (int s = 0; s < d_; ++s) { 
          // original scheme
          // S1(d_ * k + s, m) += Rcorner(nvc * s + i, nd + m);
          // S2(m, d_ * k + s) += Rcorner(nvc * s + i, nd + m);

          // mixed scheme (best one)
          // S1(d_ * k + s, m) += Rcorner(nvc * s + i, nd + m);
          // S2(m, d_ * k + s) += Ycorner(nvc * s + i, nd + m);

          // scheme based on rescaling 
          S1(d_ * k + s, m) += Rcorner(nvc * s + i, nd + m) * volume_corner[n] / u[n];
          S2(m, d_ * k + s) += Rcorner(nvc * s + i, nd + m) * volume_corner[n] / u[n];

          // support-operator
          // S1(d_ * k + s, m) += Ycorner(nvc * s + i, nd + m);
          // S2(m, d_ * k + s) += Ycorner(nvc * s + i, nd + m);
        }
      }
    }
  }

  S2 *= utot;  // only in 2D
}

}  // namespace WhetStone
}  // namespace Amanzi

