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
#include "SingleFaceMesh.hh"
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
std::vector<SchemaItem>
MFD3D_LagrangeSerendipity::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::NODE, DOF_Type::SCALAR, 1));

  if (order_ > 1) {
    int nk = PolynomialSpaceDimension(d_ - 1, order_ - 2);
    items.push_back(std::make_tuple(AmanziMesh::Entity_kind::FACE, DOF_Type::SCALAR, nk));

    if (d_ == 3) {
      nk = PolynomialSpaceDimension(d_ - 2, order_ - 2);
      items.push_back(std::make_tuple(AmanziMesh::Entity_kind::EDGE, DOF_Type::MOMENT, nk));
    }
  }

  // this should be calculated FIXME
  if (order_ > 3) {
    int nk = PolynomialSpaceDimension(d_, order_ - 4);
    items.push_back(std::make_tuple(AmanziMesh::Entity_kind::CELL, DOF_Type::SCALAR, nk));
  }

  return items;
}


/* ******************************************************************
* High-order consistency condition for the stiffness matrix.
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int
MFD3D_LagrangeSerendipity::H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, edges;
  nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  int nfaces = mesh_->getCellNumFaces(c);
  int nedges(0);
  if (d_ == 3) {
    edges = mesh_->getCellEdges(c);
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
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  // selecting regularized basis
  Basis_Regularized basis;
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

  return 0;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int
MFD3D_LagrangeSerendipity::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* Projector on a mesh face.
****************************************************************** */
void
MFD3D_LagrangeSerendipity::ProjectorFace_(int f,
                                          const std::vector<Polynomial>& ve,
                                          const ProjectorType type,
                                          const Polynomial* moments,
                                          Polynomial& uf)
{
  const auto& xf = mesh_->getFaceCentroid(f);
  const auto& normal = mesh_->getFaceNormal(f);
  SurfaceCoordinateSystem coordsys(xf, normal);

  Teuchos::RCP<SingleFaceMesh> surf_mesh =
    Teuchos::rcp(new SingleFaceMesh(mesh_, f, coordsys));
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_cache =
    Teuchos::rcp(new AmanziMesh::Mesh(surf_mesh, Teuchos::null));


  std::vector<Polynomial> vve;
  for (int i = 0; i < ve.size(); ++i) {
    Polynomial tmp(ve[i]);
    tmp.ChangeCoordinates(xf, *coordsys.tau());
    vve.push_back(tmp);
  }

  ProjectorCell_(surf_mesh_cache, 0, vve, vve, type, moments, uf);
  uf.ChangeOrigin(AmanziGeometry::Point(d_ - 1));
  uf.InverseChangeCoordinates(xf, *coordsys.tau());
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
void
MFD3D_LagrangeSerendipity::ProjectorCell_(const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
                                          int c,
                                          const std::vector<Polynomial>& ve,
                                          const std::vector<Polynomial>& vf,
                                          const ProjectorType type,
                                          const Polynomial* moments,
                                          Polynomial& uc)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->getSpaceDimension();

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mymesh, c, order_, ptmp);

  // calculate stiffness matrix
  Tensor T(d, 1);
  T(0, 0) = 1.0;

  DenseMatrix N, A;
  if (d == 2)
    H1consistency2D_(mymesh, c, T, N, A);
  else
    H1consistency3D_(c, T, N, A, false);

  // select number of non-aligned edges: we assume cell convexity
  int nfaces;
  {
    const auto& faces = mymesh->getCellFaces(c);
    nfaces = faces.size();
  }
  int eta(3);
  if (nfaces > 3) eta = 4;

  // degrees of freedom: serendipity space S contains all boundary dofs
  // plus a few internal dofs that depend on the value of eta.
  int nd = PolynomialSpaceDimension(d, order_);
  int ndof = N.NumRows();
  int ndof_c = PolynomialSpaceDimension(d, order_ - 2);
  int ndof_cs = PolynomialSpaceDimension(d, order_ - eta);
  int ndof_f(ndof - ndof_c);

  // extract submatrix
  DenseMatrix Ns, NN(nd, nd);
  Ns = N.SubMatrix(0, ndof_f, 0, nd);

  NN.Multiply(Ns, Ns, true);
  NN.InverseSPD();

  // calculate degrees of freedom (Ns^T Ns)^{-1} Ns^T v
  // for consistency with other code, we use v5 for polynomial coefficients
  const AmanziGeometry::Point& xc = mymesh->getCellCentroid(c);
  DenseVector v1(nd), v3(std::max(1, ndof_cs)), v5(nd);

  DenseVector vdof(ndof_f + ndof_cs);
  CalculateDOFsOnBoundary_(mymesh, c, ve, vf, vdof);

  // DOFs inside cell: copy moments from input data
  if (ndof_cs > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v4 = moments->coefs();
    AMANZI_ASSERT(ndof_cs == v4.NumRows());

    for (int n = 0; n < ndof_cs; ++n) { vdof(ndof_f + n) = v4(n); }
  }

  Ns.Multiply(vdof, v1, true);
  NN.Multiply(v1, v5, false);

  // this gives the least square projector
  uc = basis.CalculatePolynomial(mymesh, c, order_, v5);

  // H1 projector needs to populate moments from ndof_cs + 1 till ndof_c
  if (type == ProjectorType::H1) {
    DenseVector v4(nd);
    DenseMatrix M;
    Polynomial poly(d, order_);

    NumericalIntegration numi(mymesh);
    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

    GrammMatrix(poly, integrals_, basis, M);
    M.Multiply(v5, v4, false);

    vdof.Reshape(ndof_f + ndof_c);
    for (int n = ndof_cs; n < ndof_c; ++n) { vdof(ndof_f + n) = v4(n) / mymesh->getCellVolume(c); }

    R_.Multiply(vdof, v4, true);
    G_.Multiply(v4, v5, false);

    v5(0) = uc(0);
    uc = basis.CalculatePolynomial(mymesh, c, order_, v5);
  }

  // L2 projector is different if the set S contains some internal dofs
  if (type == ProjectorType::L2 && ndof_cs > 0) {
    DenseVector v4(nd), v6(nd - ndof_cs);
    DenseMatrix M, M2;
    Polynomial poly(d, order_);

    NumericalIntegration numi(mymesh);
    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

    GrammMatrix(poly, integrals_, basis, M);
    M2 = M.SubMatrix(ndof_cs, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    for (int n = 0; n < ndof_cs; ++n) { v4(n) = v3(n) * mymesh->getCellVolume(c); }

    for (int n = 0; n < nd - ndof_cs; ++n) { v4(ndof_cs + n) = v6(n); }

    M.InverseSPD();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mymesh, c, order_, v5);
  }

  // set correct origin
  uc.set_origin(xc);
}


/* ******************************************************************
* Calculate boundary degrees of freedom in 2D and 3D.
****************************************************************** */
void
MFD3D_LagrangeSerendipity::CalculateDOFsOnBoundary_(
  const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
  int c,
  const std::vector<Polynomial>& ve,
  const std::vector<Polynomial>& vf,
  DenseVector& vdof)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->getSpaceDimension();

  Entity_ID_List nodes;
  nodes = mymesh->getCellNodes(c);
  int nnodes = nodes.size();

  const auto& faces = mymesh->getCellFaces(c);
  int nfaces = faces.size();

  std::vector<const PolynomialBase*> polys(2);
  NumericalIntegration numi(mymesh);

  int i0, i1, pos;
  AmanziGeometry::Point xv(d);

  // number of moments of faces
  Polynomial pf;
  if (order_ > 1) { pf.Reshape(d - 1, order_ - 2); }

  int row(nnodes);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    Entity_ID_List face_nodes;
    face_nodes = mymesh->getFaceNodes(f);
    int nfnodes = face_nodes.size();

    if (d == 2) {
      for (int j = 0; j < nfnodes; j++) {
        int v = face_nodes[j];
        xv = mymesh->getNodeCoordinate(v);

        pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
        vdof(pos) = vf[n].Value(xv);
      }
    }

    if (order_ > 1) {
      double area = mymesh->getFaceArea(f);
      const AmanziGeometry::Point& xf = mymesh->getFaceCentroid(f);
      const AmanziGeometry::Point& normal = mymesh->getFaceNormal(f);

      // local coordinate system with origin at face centroid
      SurfaceCoordinateSystem coordsys(xf, normal);
      const auto& tau = *coordsys.tau();

      polys[0] = &(vf[n]);

      for (auto it = pf.begin(); it < pf.end(); ++it) {
        const int* index = it.multi_index();
        double factor = (d == 2) ? 1.0 : std::pow(area, -(double)it.MonomialSetOrder() / 2);
        Polynomial fmono(d - 1, index, factor);
        fmono.InverseChangeCoordinates(xf, tau);

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }
  }

  if (d == 3) {
    const auto& edges = mymesh->getCellEdges(c);
    int nedges = edges.size();

    Polynomial pe(d - 2, order_ - 2);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];

      // nodal DOFs
      auto nodes = mymesh->getEdgeNodes(e);
      i0 = nodes[0]; 
      i1 = nodes[1]; 

      xv = mymesh->getNodeCoordinate(i0);
      pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), i0));
      vdof(pos) = ve[n].Value(xv);

      xv = mymesh->getNodeCoordinate(i1);
      pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), i1));
      vdof(pos) = ve[n].Value(xv);

      // edge moments
      const auto& xe = mymesh->getEdgeCentroid(e);
      double length = mymesh->getEdgeLength(e);
      std::vector<AmanziGeometry::Point> tau(1, mymesh->getEdgeVector(e));

      polys[0] = &(ve[n]);

      for (auto it = pe.begin(); it < pe.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d - 2, index, 1.0);
        fmono.InverseChangeCoordinates(xe, tau);

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsEdge(e, polys) / length;
        row++;
      }
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
