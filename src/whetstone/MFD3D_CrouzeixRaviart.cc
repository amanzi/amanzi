/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_CrouzeixRaviart.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_CrouzeixRaviart::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List edges, faces, face_edges;
  std::vector<int> dirs, map;

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  N.Reshape(nedges, d_ + 1);
  Ac.Reshape(nedges, nedges);

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int num_faces = faces.size();

  // to calculate matrix R, we use temporary matrix N 
  N.PutScalar(0.0);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    if (d_ == 2) {
      for (int k = 0; k < d_; k++) N(i, k) += normal[k] * dirs[i];
    } else {
      mesh_->face_to_cell_edge_map(f, c, &map);
      int num_face_edges = map.size();

      for (int j = 0; j < num_face_edges; j++) {
        for (int k = 0; k < d_; k++) N(map[j], k) += normal[k] * dirs[i];
      }
    }
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d_), v2(d_);
  double volume = mesh_->cell_volume(c);
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; k++) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  // calculate N
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    for (int k = 0; k < d_; k++) N(i, k) = xe[k] - xc[k];
    N(i, d_) = 1;  // additional column is added to the consistency condition
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. 
****************************************************************** */
int MFD3D_CrouzeixRaviart::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* High-order consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_CrouzeixRaviart::H1consistencyHO(
    int c, int order, const Tensor& K,
    DenseMatrix& N, DenseMatrix& R, DenseMatrix& Ac, DenseMatrix& G)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c); 

  // calculate degrees of freedom 
  Polynomial poly(d_, order), pf(d_ - 1, order - 1), pc;
  if (order > 1) {
    pc.Reshape(d_, order - 2);
  }
  int nd = poly.size();
  int ndf = pf.size();
  int ndc = pc.size();

  int ndof = nfaces * ndf + ndc;
  N.Reshape(ndof, nd);
  R.Reshape(ndof, nd);
  Ac.Reshape(ndof, ndof);
  G.Reshape(nd, nd);

  // pre-calculate integrals of monomials 
  NumericalIntegration numi(mesh_);
  Polynomial integrals(d_, 2 * order - 2);

  for (int k = 0; k <= 2 * order - 2; ++k) {
    numi.IntegrateMonomialsCell(c, integrals.monomials(k));
  }

  // populate matrices N and R
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  R.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const Polynomial*> polys(2);

  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) { 
    const int* index = it.multi_index();
    Polynomial cmono(d_, index);
    cmono.set_origin(xc);  

    // N and R: degrees of freedom on faces 
    VectorPolynomial grad;
    grad.Gradient(cmono);
     
    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(0);

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[i];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 

      // local coordinate system with origin at face centroid
      if (d_ == 2) {
        tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
      }

      Polynomial tmp = grad * normal;
      tmp.ChangeCoordinates(xf, tau);

      for (auto jt = tmp.begin(); jt.end() <= tmp.end(); ++jt) {
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        int n = jt.PolynomialPosition();
        R(row + n, col) = tmp(m, k);
      }

      for (auto jt = pf.begin(); jt.end() <= pf.end(); ++jt) {
        const int* jndex = jt.multi_index();
        Polynomial fmono(d_ - 1, jndex);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        int n = jt.PolynomialPosition();
        N(row + n, col) = numi.IntegratePolynomialsFace(f, polys);
      }
      row += ndf;
    }

    // N and R: degrees of freedom in cells
    if (cmono.order() > 1) {
      Polynomial tmp = cmono.Laplacian();
      for (auto jt = tmp.begin(); jt.end() <= tmp.end(); ++jt) {
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        int n = jt.PolynomialPosition();

        R(row + n, col) = -tmp(m, k);

        int nm(0);
        const int* jndex = jt.multi_index();
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        double sum(0.0), tmp;
        for (int i = 0; i < d_; ++i) {
          if (index[i] > 1) {
            multi_index[i] -= 2;
            const auto& coefs = integrals.monomials(nm - 2).coefs();
            tmp = coefs[poly.MonomialPosition(multi_index)]; 
            sum += tmp * index[i] * (index[i] - 1);
            multi_index[i] += 2;
          }
        }

        N(row + n, col) = K(0, 0) * sum;
      }
    }
  }

  // set the Gramm-Schidt matrix for polynomials
  G.PutScalar(0.0);

  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt.end() <= poly.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      
      int n(0);
      int multi_index[3];
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = index[i] + jndex[i];
        n += multi_index[i];
      }

      double sum(0.0), tmp;
      for (int i = 0; i < d_; ++i) {
        if (index[i] > 0 && jndex[i] > 0) {
          multi_index[i] -= 2;
          const auto& coefs = integrals.monomials(n - 2).coefs();
          tmp = coefs[poly.MonomialPosition(multi_index)]; 
          sum += tmp * index[i] * jndex[i];
          multi_index[i] += 2;
        }
      }

      G(l, k) = G(k, l) = K(0, 0) * sum; 
    }
  }

  // calculate R inv(G) R^T
  DenseMatrix RG(ndof, nd), Rtmp(nd, ndof);

  G.Inverse();
  RG.Multiply(R, G, false);

  Rtmp.Transpose(R);
  Ac.Multiply(RG, Rtmp, false);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_CrouzeixRaviart::StiffnessMatrixHO(
    int c, int order, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N, R, G;

  int ok = H1consistencyHO(c, order, K, N, R, A, G);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



