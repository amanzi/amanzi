/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Serendipity CrouzeixRaviart-type element: degrees of freedom are 
  moments on faces and inside cell. The number of later is 
  reduced significantly for polytopal cells. 
*/

#include <cmath>
#include <tuple>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "Basis_Regularized.hh"
#include "SurfaceCoordinateSystem.hh"
#include "GrammMatrix.hh"
#include "MFD3D_CrouzeixRaviartSerendipity.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> MFD3D_CrouzeixRaviartSerendipity::schema() const
{
  int nk = PolynomialSpaceDimension(d_ - 1, order_ - 1);
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::MOMENT, nk));

  if (order_ > 3) {
    nk = PolynomialSpaceDimension(d_, order_ - 4);
    items.push_back(std::make_tuple(AmanziMesh::CELL, DOF_Type::MOMENT, nk));
  }

  return items;
}


/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_CrouzeixRaviartSerendipity::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  int nfaces = mesh_->cell_get_num_faces(c);
  AMANZI_ASSERT(nfaces > 3);  // FIXME

  // calculate degrees of freedom 
  Polynomial poly(d_, order_), pf, pc;
  pf.Reshape(d_ - 1, order_ - 1);
  if (order_ > 3) {
    pc.Reshape(d_, order_ - 4);
  }

  int nd = poly.size();
  int ndf = pf.size();
  int ndc = pc.size();
  int ndof_S = nfaces * ndf + ndc;

  // calculate full matrices
  DenseMatrix Nf, Af;
  MFD3D_CrouzeixRaviartAnyOrder::H1consistency(c, K, Nf, Af);

  // pre-calculate integrals of monomials 
  NumericalIntegration<AmanziMesh::StaticMesh> numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  // selecting regularized basis
  Basis_Regularized<AmanziMesh::StaticMesh> basis;
  basis.Init(mesh_, c, order_, integrals_.poly());

  // Gramm matrix for polynomials
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

  return 0;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_CrouzeixRaviartSerendipity::StiffnessMatrix(
    int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* L2 projector 
****************************************************************** */
void MFD3D_CrouzeixRaviartSerendipity::ProjectorCell_(
    int c, const std::vector<Polynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, Polynomial& uc)
{
  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized<AmanziMesh::StaticMesh> basis;
  basis.Init(mesh_, c, order_, ptmp);

  // calculate full matrices
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  MFD3D_CrouzeixRaviartAnyOrder::H1consistency(c, T, N, A);  

  // select number of non-aligned edges: we assume cell convexity 
  int nfaces = mesh_->cell_get_num_faces(c);
  int eta(3);
  if (nfaces > 3) eta = 4;

  // number of degrees of freedom
  int nd = G_.NumRows();
  int ndof = A.NumRows();
  int ndof_c = PolynomialSpaceDimension(d_, order_ - 2);
  int ndof_cs = PolynomialSpaceDimension(d_, order_ - eta);
  int ndof_f(ndof - ndof_c);

  // extract submatrix
  DenseMatrix Ns, NN(nd, nd);
  Ns = N.SubMatrix(0, ndof_f, 0, nd);

  NN.Multiply(Ns, Ns, true);
  NN.Inverse();

  // calculate degrees of freedom
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  DenseVector v1(nd), v5(nd);

  DenseVector vdof(ndof_f + ndof_cs);
  CalculateDOFsOnBoundary_(c, vf, vdof);

  // DOFs inside cell: copy moments from input data
  if (ndof_cs > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    for (int n = 0; n < ndof_cs; ++n) {
      vdof(ndof_f + n) = v3(n);
    }
  }

  Ns.Multiply(vdof, v1, true);
  NN.Multiply(v1, v5, false);

  uc = basis.CalculatePolynomial(mesh_, c, order_, v5);

  // H1 projector needs to populate moments from ndof_cs + 1 till ndof_c
  if (type == ProjectorType::H1) {
    DenseVector v4(nd);
    DenseMatrix M;
    Polynomial poly(d_, order_);

    NumericalIntegration<AmanziMesh::StaticMesh> numi(mesh_);
    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

    GrammMatrix(poly, integrals_, basis, M);
    M.Multiply(v5, v4, false);

    vdof.Reshape(ndof_f + ndof_c);
    for (int n = ndof_cs; n < ndof_c; ++n) {
      vdof(ndof_f + n) = v4(n) / mesh_->cell_volume(c); 
    }

    R_.Multiply(vdof, v4, true);
    G_.Multiply(v4, v5, false);

    v5(0) = uc(0);
    uc = basis.CalculatePolynomial(mesh_, c, order_, v5);
  }

  // L2 projector is different if the set S contains some internal dofs
  if (type == ProjectorType::L2 && ndof_cs > 0) {
    DenseVector v4(nd), v6(nd - ndof_cs);
    DenseMatrix M, M2;
    Polynomial poly(d_, order_);

    NumericalIntegration<AmanziMesh::StaticMesh> numi(mesh_);
    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

    GrammMatrix(poly, integrals_, basis, M);
    M2 = M.SubMatrix(ndof_cs, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    const DenseVector& v3 = moments->coefs();
    for (int n = 0; n < ndof_cs; ++n) {
      v4(n) = v3(n) * mesh_->cell_volume(c);
    }

    for (int n = 0; n < nd - ndof_cs; ++n) {
      v4(ndof_cs + n) = v6(n);
    }

    M.Inverse();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mesh_, c, order_, v5);
  }

  // set correct origin
  uc.set_origin(xc);
}


/* ******************************************************************
* Calculate degrees of freedom in 2D.
****************************************************************** */
void MFD3D_CrouzeixRaviartSerendipity::CalculateDOFsOnBoundary_(
    int c, const std::vector<Polynomial>& vf, DenseVector& vdof)
{
  Entity_ID_List faces;

  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  std::vector<const PolynomialBase*> polys(2);
  NumericalIntegration<AmanziMesh::StaticMesh> numi(mesh_);

  AmanziGeometry::Point xv(d_);

  // number of moments of faces
  Polynomial pf(d_ - 1, order_ - 1);

  int row(0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    // local coordinate system with origin at face centroid
    SurfaceCoordinateSystem coordsys(xf, normal);

    polys[0] = &(vf[n]);

    for (auto it = pf.begin(); it < pf.end(); ++it) {
      const int* index = it.multi_index();
      Polynomial fmono(d_ - 1, index, 1.0);
      fmono.InverseChangeCoordinates(xf, *coordsys.tau());  

      polys[1] = &fmono;

      vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
      row++;
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

