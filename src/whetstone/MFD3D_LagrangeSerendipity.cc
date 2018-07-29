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

#include "Basis_Regularized.hh"
#include "CoordinateSystems.hh"
#include "GrammMatrix.hh"
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

  // selecting regularized basis
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_);

  // Dot-product matrix for polynomials and Laplacian of polynomials
  DenseMatrix M(nd, nd);
  GrammMatrix(poly, integrals_, basis, M);

  // setup matrix representing Laplacian of polynomials
  DenseMatrix L(nd, nd);
  L.PutScalar(0.0);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d_, index, factor);
    VectorPolynomial grad, Kgrad;

    grad.Gradient(cmono);
    Kgrad = K * grad;
    Polynomial lap = Divergence(Kgrad);
    
    for (auto jt = lap.begin(); jt < lap.end(); ++jt) {
      int l = jt.PolynomialPosition();
      int m = jt.MonomialSetOrder();
      int n = jt.MonomialSetPosition();
      L(l, k) = lap(m, n) / basis.monomial_scales()[m];
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


/* ******************************************************************
* L2 projector 
****************************************************************** */
void MFD3D_LagrangeSerendipity::L2Cell(
    int c, const std::vector<VectorPolynomial>& vf,
    VectorPolynomial& moments, VectorPolynomial& uc)
{
  // create integration object for a single cell
  NumericalIntegration numi(mesh_);

  // selecting regularized basis
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  MFD3D_Lagrange::H1consistency(c, T, N, A);  

  // number of degrees of freedom
  Polynomial pc;
  if (order_ > 1) {
    pc.Reshape(d_, order_ - 2);
  }

  int nd = G_.NumRows();
  int ndof = A.NumRows();
  int ndof_c(pc.size());
  int ndof_f(ndof - ndof_c);

  // extract submatrix
  DenseMatrix Ns, NN(nd, nd);
  Ns = N.SubMatrix(0, ndof_f, 0, nd);

  NN.Multiply(Ns, Ns, true);
  NN.Inverse();

  // calculate degrees of freedom (Ns^T Ns)^{-1} Ns^T v
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  DenseVector vdof(ndof_f), v1(nd), v2(nd);

  int dim = vf[0].size();
  uc.resize(dim);

  for (int i = 0; i < dim; ++i) {
    CalculateDOFsOnBoundary_(c, vf, vdof, i);

    Ns.Multiply(vdof, v1, true);
    NN.Multiply(v1, v2, false);

    uc[i] = basis.CalculatePolynomial(mesh_, c, order_, v2);

    // set correct origin 
    uc[i].set_origin(xc);
  }
}


/* ******************************************************************
* Calculate degrees of freedom in 2D.
****************************************************************** */
void MFD3D_LagrangeSerendipity::CalculateDOFsOnBoundary_(
    int c, const std::vector<VectorPolynomial>& vf,
    DenseVector& vdof, int i)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  std::vector<const Polynomial*> polys(2);
  NumericalIntegration numi(mesh_);

  AmanziGeometry::Point xv(d_);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);

  // number of moments of faces
  Polynomial pf;
  if (order_ > 1) {
    pf.Reshape(d_ - 1, order_ - 2);
  }

  int row(nnodes);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int nfnodes = face_nodes.size();

    for (int j = 0; j < nfnodes; j++) {
      int v = face_nodes[j];
      mesh_->node_get_coordinates(v, &xv);

      int pos = WhetStone::FindPosition(v, nodes);
      vdof(pos) = vf[n][i].Value(xv);
    }

    if (order_ > 1) { 
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      FaceCoordinateSystem(normal, tau);

      polys[0] = &(vf[n][i]);

      for (auto it = pf.begin(); it < pf.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

