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
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

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
  DenseMatrix Nf, Rf, Gf, Af;
  MFD3D_Lagrange::H1consistencyHO(c, order_, K, Nf, Rf, Gf, Af);

  // calculate submatrix
  DenseMatrix R, M;
  N = Nf.SubMatrix(0, ndof_S, 0, nd);
  R.Reshape(ndof_S, nd);
  Ac.Reshape(ndof_S, ndof_S);

  // pre-calculate integrals of monomials 
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  // set the Gramm-Schidt matrix for polynomials
  M.Reshape(nd, nd);
  M.PutScalar(0.0);

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

  // calculate R inv(G) R^T
  DenseMatrix RG(ndof_S, nd), Rtmp(nd, ndof_S);

  RG.Multiply(R, Gf, false);
  Rtmp.Transpose(R);
  Ac.Multiply(RG, Rtmp, false);

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



