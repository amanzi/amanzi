/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Serendipity Lagrange-type element: degrees of freedom are nodal values,
  moments on edges, and selected moments on faces and inside cell:
  Degrees of freedom are ordered as follows:
    (1) nodal values in the natural order;
    (2) moments on faces groupped by face;
    (3) moments on edges, groupped by edge (in 3D);
    (4) moments inside cell.
*/

#include <cmath>
#include <tuple>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "Basis_Regularized.hh"
#include "GrammMatrix.hh"
#include "MFD3D_LagrangeAnyOrder.hh"
#include "MFD3D_LagrangeSerendipity.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
MFD3D_LagrangeSerendipity::MFD3D_LagrangeSerendipity(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D_LagrangeAnyOrder(plist, mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> MFD3D_LagrangeSerendipity::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::NODE, DOF_Type::SCALAR, 1));

  if (order_ > 1) {
    int nk = PolynomialSpaceDimension(d_ - 1, order_ - 2);
    items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, nk));
  }

  // this should be calculated FIXME
  if (order_ > 3) {
    int nk = PolynomialSpaceDimension(d_, order_ - 4);
    items.push_back(std::make_tuple(AmanziMesh::CELL, DOF_Type::SCALAR, nk));
  }

  return items;
}


/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_LagrangeSerendipity::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, edges;
  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  int nfaces = mesh_->cell_get_num_faces(c);
  int nedges(0);
  if (d_ == 3) {
    mesh_->cell_get_edges(c, &edges);
    nedges = edges.size();
  }

  // select number of non-aligned edges: we assume cell convexity 
  int eta(3);
  if (nfaces > 3) eta = 4;

  // calculate degrees of freedom: serendipity space S contains all boundary
  // dofs plus a few internal dofs that depedend on the value of eta.
  Polynomial poly(d_, order_);

  int nd, nde(0), ndf, ndc;
  nd = poly.size();
  nde = (d_ == 2) ? 0 : PolynomialSpaceDimension(d_ - 2, order_ - 2);
  ndf = PolynomialSpaceDimension(d_ - 1, order_ - 2);
  ndc = PolynomialSpaceDimension(d_, order_ - eta);
  int ndof_S = nnodes + nedges * nde + nfaces * ndf + ndc;

  // calculate full matrices
  DenseMatrix Nf, Af;
  MFD3D_LagrangeAnyOrder::H1consistency(c, K, Nf, Af);

  // pre-calculate integrals of monomials 
  NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  // selecting regularized basis
  Basis_Regularized<AmanziMesh::Mesh> basis;
  basis.Init(mesh_, c, order_, integrals_.poly());

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

    auto Kgrad = K * Gradient(cmono);
    Polynomial lap = Divergence(Kgrad);
    
    for (auto jt = lap.begin(); jt < lap.end(); ++jt) {
      int l = jt.PolynomialPosition();
      int m = jt.MonomialSetOrder();
      L(l, k) = lap(l) / basis.monomial_scales()[m];
    }  
  }

  // calculate matrices N and R
  // -- extract sub-matrices. More work is required in 3D
  DenseMatrix Rf(R_);
  N = Nf.SubMatrix(0, ndof_S, 0, nd);
  R_ = Rf.SubMatrix(0, ndof_S, 0, nd);

  // -- add correcton Ns (Ns^T Ns)^{-1} M L to matrix R_
  DenseMatrix NN(nd, nd), NM(nd, nd);

  NN.Multiply(N, N, true);
  NN.InverseSPD();

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
* Calculate boundary degrees of freedom in 2D and 3D.
****************************************************************** */
void MFD3D_LagrangeSerendipity::CalculateDOFsOnBoundary_(
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf, DenseVector& vdof)
{
  Entity_ID_List nodes, faces, edges;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  std::vector<const PolynomialBase*> polys(2);
  NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);

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

      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
      vdof(pos) = vf[n].Value(xv);
    }

    if (order_ > 1) { 
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      SurfaceCoordinateSystem coordsys(normal);
      const auto& tau = *coordsys.tau();

      polys[0] = &(vf[n]);

      for (auto it = pf.begin(); it < pf.end(); ++it) {
        const int* index = it.multi_index();
        double factor = (d_ == 2) ? 1.0 : std::pow(area, -(double)it.MonomialSetOrder() / 2);
        Polynomial fmono(d_ - 1, index, factor);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }
  }

  if (d_ == 3) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    Polynomial pe(d_ - 2, order_ - 2);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      const auto& xe = mesh_->edge_centroid(e);
      double length = mesh_->edge_length(e);
      std::vector<AmanziGeometry::Point> tau(1, mesh_->edge_vector(e));

      polys[0] = &(ve[n]);

      for (auto it = pe.begin(); it < pe.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 2, index, 1.0);
        fmono.InverseChangeCoordinates(xe, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsEdge(e, polys) / length;
        row++;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

