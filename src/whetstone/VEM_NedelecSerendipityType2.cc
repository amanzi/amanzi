/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  High-order 3D serendipity Nedelec type 2 element: degrees of freedom
  are moments on edges, selected moments on faces and inside cell:
  Degrees of freedom are ordered as follows:
    (1) moments on edges, order is moment number -> edge id
    (2) moments on faces, order is moment number -> face id;
    (4) moments inside cell
  Vector degrees of freedom are ordered first by moments on a geometric
  entity and then by the vector component.

  At the moment, loop over the space of test polynomials is hard-coded.
*/

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "Basis_Regularized.hh"
#include "GrammMatrix.hh"
#include "VEM_NedelecSerendipityType2.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
VEM_NedelecSerendipityType2::VEM_NedelecSerendipityType2(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D(mesh),
    DeRham_Edge(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> VEM_NedelecSerendipityType2::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::EDGE, DOF_Type::SCALAR, order_));

  return items;
}


/* ******************************************************************
* VEM scheme:: mass matrix for second-order scheme.
****************************************************************** */
int VEM_NedelecSerendipityType2::L2consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  AMANZI_ASSERT(d_ == 3);

  Entity_ID_List edges, faces, fedges;
  std::vector<int> fdirs, edirs;

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  double volume = mesh_->cell_volume(c);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  // selecting regularized basis (parameter integrals_ is not used)
  Basis_Regularized<AmanziMesh::Mesh> basis;
  basis.Init(mesh_, c, order_, integrals_.poly());

  // pre-calculate integrals of monomials 
  NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  // calculate degrees of freedom: serendipity space S contains all boundary dofs
  Polynomial pc(d_, order_), pf(d_ - 1, order_), pe(d_ - 2, order_);

  int ndc, ndf, nde;
  ndc = pc.size();
  ndf = PolynomialSpaceDimension(d_ - 1, order_);
  nde = PolynomialSpaceDimension(d_ - 2, order_);
  int ndof_S = nedges * nde;

  // iterators
  VectorPolynomialIterator it0(d_, d_, order_), it1(d_, d_, order_);
  it0.begin();
  it1.end();

  // H-curl spaces on faces
  std::vector<WhetStone::DenseMatrix> vNf, vMG;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
    auto surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));

    DenseMatrix Nf, Mf, MG;
    L2consistency2D_<SurfaceMiniMesh>(surf_mesh, f, K, Nf, Mf, MG);

    vNf.push_back(Nf);
    vMG.push_back(MG);
    vcoordsys.push_back(coordsys);
  }

  // contributions of least square projectors on faces
  Tensor Id(d_ - 1, 2);
  Id.MakeDiagonal(1.0);

  N.Reshape(ndof_S, ndc * d_);
  N.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    const auto& Nf = vNf[n];
    DenseMatrix NfT(Nf.NumCols(), Nf.NumRows());
    NfT.Transpose(Nf);

    auto NN = NfT * Nf;
    NN.InverseMoorePenrose();

    // matrix Nf * NN represents the least square projector
    auto P0 = Nf * NN * (vMG[n] ^ Id);

    mesh_->face_get_edges_and_dirs(f, &fedges, &edirs);
    int nfedges = fedges.size();

    // local face -> local cell map
    std::vector<int> map;
    for (int i = 0; i < nfedges; ++i) {
      int e = fedges[i];
      int pos = std::distance(edges.begin(), std::find(edges.begin(), edges.end(), e));
      for (int l = 0; l < nde; ++l) map.push_back(nde * pos + l); 
    }

    // integral over face requires vector of polynomial coefficients of pc
    WhetStone::DenseVector v(P0.NumCols()), p0v(P0.NumRows());
    auto& tau = *vcoordsys[n]->tau();

    for (auto it = it0; it < it1; ++it) {
      int k = it.VectorComponent();
      int m = it.PolynomialPosition();

      const int* index = it.multi_index();
      double factor = basis.monomial_scales()[it.MonomialSetOrder()];
      Polynomial qe(d_, index, factor);
      qe.set_origin(xc);
      qe.ChangeCoordinates(xf, tau);

      int s(0);
      /*
      for (int i = 0; i < d_ - 1; ++i) {
        Polynomial q1(d_, 1);
        for (int l = 0; l < d_; ++l) {
          q1(l + 1) = tau[i][l] * normal[k] - tau[i][k] * normal[l];
        }
        q1.set_origin(xc);
        q1.ChangeCoordinates(xf, tau);

        q1 *= qe;
      */

      AmanziGeometry::Point p(d_);
      p[k] = 1.0 / area;
      auto tmp = vcoordsys[n]->Project(((xf - xc) ^ p) ^ normal, false);

      for (int l = 0; l < d_ - 1; ++l) v(l) = tmp[l];
      P0.Multiply(v, p0v, false);

      for (int i = 0; i < p0v.NumRows(); ++i) {
        N(map[i], d_ * m + k) += p0v(i) * fdirs[n] / 2;
      }
    }
  }

  // calculate Mc = R (R^T N)^{-1} R^T 
  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);

  Tensor Kinv(K);
  Kinv.Inverse();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = Kinv * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = N(j, k);
      Mc(i, j) = (v2 * v3) / volume;
    }
  }

  // Rows of matrix N are simply tangents. Since N goes to the Gramm-Schmidt 
  // orthogonalizetion procedure, we drop scaling with tensorial factor K.
  int row(0);
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);

    for (int k = 0; k < d_; ++k) N(i, k) = tau[k] / len;
  }

  // create Gramm-Schmidt for serial and vector cases following the order of dofs
  /*
  DenseMatrix M1(nd, nd);
  GrammMatrix(poly, integrals_, basis, M1);

  Tensor Kinv(d_, 2);
  if (K.rank() == 2) Kinv = K;
  else Kinv += K(0, 0);
  Kinv.Inverse();

  auto M3 = M1 ^ Kinv;

  // Integrate least-square polynomial
  DenseMatrix P0T(P0.NumCols(), P0.NumRows());
  P0T.Transpose(P0);
  Mc = P0 * M3 * P0T;
  */

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix for edge-based discretization.
****************************************************************** */
int VEM_NedelecSerendipityType2::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, K, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi

