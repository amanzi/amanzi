/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Lagrange-type element: degrees of freedom are nodal values and
  moments on edges, faces and inside cell.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Lagrange.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_Lagrange::H1consistencyHO(
    int c, int order, const Tensor& K,
    DenseMatrix& N, DenseMatrix& R, DenseMatrix& G, DenseMatrix& Ac)
{
  ASSERT(d_ == 2);  // FIXME
  ASSERT(order < 3);

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c); 
  double volume = mesh_->cell_volume(c); 

  // calculate degrees of freedom 
  Polynomial poly(d_, order), pf, pc;
  if (order > 1) {
    pf.Reshape(d_ - 1, order - 2);
    pc.Reshape(d_, order - 2);
  }
  int nd = poly.size();
  int ndf = pf.size();
  int ndc = pc.size();

  int ndof = nnodes + nfaces * ndf + ndc;
  N.Reshape(ndof, nd);
  R.Reshape(ndof, nd);
  Ac.Reshape(ndof, ndof);
  G.Reshape(nd, nd);

  // pre-calculate integrals of monomials 
  NumericalIntegration numi(mesh_);
  integrals_.Reshape(d_, 2 * order - 2, true);

  for (int k = 0; k <= 2 * order - 2; ++k) {
    numi.IntegrateMonomialsCell(c, integrals_.monomials(k));
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
    int row(nnodes);

    AmanziGeometry::Point xv(d_);
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      mesh_->node_get_coordinates(v, &xv);
      N(i, col) = cmono.Value(xv);
    }

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      AmanziGeometry::Point normal = mesh_->face_normal(f);
      if (d_ == 2) {
        tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
      }
      normal *= dirs[i];

      Entity_ID_List face_nodes;
      mesh_->face_get_nodes(f, &face_nodes);
      int nfnodes = face_nodes.size();

      if (order == 1 && col > 0) {
        for (int j = 0; j < nfnodes; j++) {
          int v = face_nodes[j];
          int pos = FindPosition(v, nodes);
          R(pos, col) += normal[col - 1] * dirs[i] / 2;
        }
      } else if (order == 2) {
        Polynomial tmp = grad * normal;

        // Simpson rule with 3 points
        int v = face_nodes[0];
        int pos0 = FindPosition(v, nodes);
        mesh_->node_get_coordinates(v, &xv);
        double q0 = tmp.Value(xv);

        v = face_nodes[1];
        int pos1 = FindPosition(v, nodes);
        mesh_->node_get_coordinates(v, &xv);
        double q1 = tmp.Value(xv);

        double qmid = tmp.Value(mesh_->face_centroid(f));

        R(pos0, col) += (q0 - qmid) / 6;
        R(pos1, col) += (q1 - qmid) / 6;
        R(row,  col) = qmid;

        for (auto jt = pf.begin(); jt.end() <= pf.end(); ++jt) {
          const int* jndex = jt.multi_index();
          Polynomial fmono(d_ - 1, jndex);
          fmono.InverseChangeCoordinates(xf, tau);  

          polys[1] = &fmono;

          int n = jt.PolynomialPosition();
          N(row + n, col) = numi.IntegratePolynomialsFace(f, polys) / area;
        }
        row += ndf;
      }
    }

    // N and R: degrees of freedom in cells
    if (cmono.order() > 1) {
      Polynomial tmp = cmono.Laplacian();
      for (auto jt = tmp.begin(); jt.end() <= tmp.end(); ++jt) {
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        int n = jt.PolynomialPosition();

        R(row + n, col) = -tmp(m, k) * volume;
      }
    }

    if (order > 1) {
      for (auto jt = pc.begin(); jt.end() <= pc.end(); ++jt) {
        int n = jt.PolynomialPosition();
        const int* jndex = jt.multi_index();

        int nm(0);
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        const auto& coefs = integrals_.monomials(nm).coefs();
        N(row + n, col) = coefs[poly.MonomialPosition(multi_index)] / volume; 
      }
    }
  }

  // set the Gramm-Schidt matrix for gradients of polynomials
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
          const auto& coefs = integrals_.monomials(n - 2).coefs();
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

  // to invert generate matrix, we add and subtruct positive number
  G(0, 0) = 1.0;
  G.Inverse();
  G(0, 0) = 0.0;
  RG.Multiply(R, G, false);

  Rtmp.Transpose(R);
  Ac.Multiply(RG, Rtmp, false);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_Lagrange::StiffnessMatrixHO(
    int c, int order, const Tensor& K,
    DenseMatrix& R, DenseMatrix& G, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistencyHO(c, order, K, N, R, G, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



