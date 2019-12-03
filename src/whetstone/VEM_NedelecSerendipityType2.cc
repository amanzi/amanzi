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
#include "VectorObjectsUtils.hh"
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

  std::vector<const PolynomialBase*> polys(2);

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

  // fixed vector
  VectorPolynomial xyz(d_, d_, 1);
  for (int i = 0; i < d_; ++i) xyz[i](i + 1) = 1.0;
  xyz.set_origin(xc);

  // Rows of matrix N are simply tangents. Since N goes to the Gramm-Schmidt 
  // orthogonalizetion procedure, we drop scaling with tensorial factor K.
  N.Reshape(ndof_S, ndc * d_);

  for (auto it = it0; it < it1; ++it) { 
    int k = it.VectorComponent();
    int col = it.VectorPolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d_, index, factor);
    cmono.set_origin(xc);  

    int row(0);
    for (int i = 0; i < nedges; i++) {
      int e = edges[i];
      const auto& xe = mesh_->edge_centroid(e);
      const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
      double length = mesh_->edge_length(e);
      std::vector<AmanziGeometry::Point> tau_edge(1, tau);

      // order 0 scheme
      // for (int k = 0; k < d_; ++k) N(i, k) = tau[k] / length;

      for (auto it = pe.begin(); it < pe.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial emono(d_ - 2, index, tau[k]);
        emono.InverseChangeCoordinates(xe, tau_edge);  

        polys[0] = &cmono;
        polys[1] = &emono;

        N(row, col) = numi.IntegratePolynomialsEdge(e, polys) / length;
        row++;
      }
    }
  }

  // ----------------------
  // pre-compute projectors
  // ----------------------
  // -- L2 projectors on faces
  Polynomial qf(d_ - 1, order_ + 1);
  std::vector<WhetStone::DenseMatrix> vL2f, vMGf;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

  Tensor Idf(d_ - 1, 2), Idc(d_, 2);
  Idf.MakeDiagonal(1.0);
  Idc.MakeDiagonal(1.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
    auto surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));

    DenseMatrix Nf, NfT, Mf, MG;
    L2consistency2D_<SurfaceMiniMesh>(surf_mesh, f, K, Nf, Mf, MG);

    vcoordsys.push_back(coordsys);

    NfT.Transpose(Nf);
    auto NN = NfT * Nf;
    NN.InverseMoorePenrose();
    auto L2f = (MG ^ Idf) * NN * NfT;
    vL2f.push_back(L2f);

    // GrammMatrix(qf, integrals_, basis, MG);
  }

  // -- L2-projector in cell
  DenseMatrix MG(ndc, ndc), NT, L2c;
  GrammMatrix(numi, order_, integrals_, basis, MG);

  NT.Transpose(N);
  auto NN = NT * N;
  NN.InverseMoorePenrose();
  L2c = (MG ^ Idc) * NN * NT;

  // contributions of least square projectors on faces
  DenseMatrix R(ndof_S, ndc * d_);
  R.PutScalar(0.0);

  for (auto it = it0; it < it1; ++it) {
    int k = it.VectorComponent();
    int m = it.PolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Monomial q(d_, index, factor);
    q.set_origin(xc);

    Polynomial cmono(q);

    // vector decomposition of vector (0, q, 0) with q at k-th position
    VectorPolynomial p1;
    Polynomial p2;
    VectorDecomposition3DCurl(q, k, p1, p2);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

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
      WhetStone::DenseVector p0v(vL2f[n].NumCols());
      auto& tau = *vcoordsys[n]->tau();
      cmono.ChangeCoordinates(xf, tau);

      VectorPolynomial p3D = p1 * (xyz * normal) - xyz * (p1 * normal);
      VectorPolynomial p2D = ProjectVectorPolynomialOnManifold(p3D, xf, *vcoordsys[n]->tau());
      auto v = ExpandCoefficients(p2D);
      v /= area;

      vL2f[n].Multiply(v, p0v, true);

      for (int i = 0; i < p0v.NumRows(); ++i) {
        R(map[i], d_ * m + k) += p0v(i) * fdirs[n];
      }
    }
  }

  // calculate Mc = R (R^T N)^{-1} R^T 
  AmanziGeometry::Point v1(d_), v2(d_), v3(d_);

  Tensor Kinv(K);
  Kinv.Inverse();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = R(i, k);
    v2 = Kinv * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = R(j, k);
      Mc(i, j) = (v2 * v3) / volume;
    }
  }

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

