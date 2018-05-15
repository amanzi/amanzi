/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Serendipity Lagrange-type element: degrees of freedom are nodal values 
  and moments on edges, faces and inside cell. The number of later is 
  reduced significantly for polygonal cells. 
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_LagrangeSerendipity::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  int nfaces = mesh_->cell_get_num_faces(c);
  AMANZI_ASSERT(nfaces > 3);  // FIXME

  // calculate degrees of freedom 
  Polynomial poly(d_, order_), pf, pc;
  if (order_ > 1) {
    pf.Reshape(d_ - 1, order_ - 2);
  }
  if (order_ > 3) {
    pc.Reshape(d_, order_ - 4);
  }

  int nd = poly.size();
  int ndf = pf.size();
  int ndc = pc.size();
  int ndof_S = nnodes + nfaces * ndf + ndc;

  // calculate full matrices
  DenseMatrix Nf, Af;
  MFD3D_Lagrange::H1consistency(c, K, Nf, Af);

  // pre-calculate integrals of monomials 
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  // set the Gramm-Schidt matrix for polynomials
  DenseMatrix M(nd, nd);

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
      const auto& coefs = integrals_.poly().monomials(n).coefs();
      M(l, k) = M(k, l) = coefs[poly.MonomialPosition(multi_index)]; 
    }
  }

  // setup matrix representing Laplacian of polynomials
  double volume = mesh_->cell_volume(c); 
  double scale = numi.MonomialNaturalScale(1, volume);

  DenseMatrix L(nd, nd);
  L.PutScalar(0.0);

  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    double factor = numi.MonomialNaturalScale(it.MonomialOrder(), volume);
    Polynomial mono(d_, index, factor);
    Polynomial lap = mono.Laplacian();
    
    for (auto jt = lap.begin(); jt.end() <= lap.end(); ++jt) {
      int l = jt.PolynomialPosition();
      int m = jt.MonomialOrder();
      int n = jt.MonomialPosition();
      L(l, k) = lap(m, n) / std::pow(scale, std::min(2, m));
    }  
  }

  // calculate matrices N and R
  // -- extract sub-matrices 
  DenseMatrix Rf(R_);
  N = Nf.SubMatrix(0, ndof_S, 0, nd);
  R_ = Rf.SubMatrix(0, ndof_S, 0, nd);

  // -- add correcton Ns (Ns^T Ns)^{-1} M L to matrix R_
  DenseMatrix NN(nd, nd), NM(nd, nd);

  NN.Multiply(N, N, true);
  NN.Inverse();

  NM.Multiply(NN, M, false);
  NN.Multiply(NM, L, false);

  Nf.Reshape(ndof_S, nd);
  Nf.Multiply(N, NN, false);

  R_ -= Nf;

  // calculate Ac = R inv(G) R^T
  Ac.Reshape(ndof_S, ndof_S);
  DenseMatrix Rtmp(nd, ndof_S);

  Nf.Multiply(R_, G_, false);
  Rtmp.Transpose(R_);
  Ac.Multiply(Nf, Rtmp, false);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_LagrangeSerendipity::StiffnessMatrix(
    int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



