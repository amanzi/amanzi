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
#include "FunctionPower.hh"
#include "GrammMatrix.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "SurfaceMiniMesh.hh"
#include "Tensor.hh"
#include "VectorObjectsUtils.hh"
#include "VEM_NedelecSerendipityType2.hh"

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
  if (d_ == 2) {
    DenseMatrix MG;
    L2consistency2D_<AmanziMesh::Mesh>(mesh_, c, K, N, Mc, MG);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }

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
  basis.Init(mesh_, c, order_ + 1, integrals_.poly());

  // pre-calculate integrals of monomials 
  NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  std::vector<const WhetStoneFunction*> funcs(2);

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
    int m = it.MonomialSetOrder();
    int k = it.VectorComponent();
    int col = it.VectorPolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[m];
    Monomial cmono(d_, index, factor);
    cmono.set_origin(xc);  

    int row(0);
    for (int i = 0; i < nedges; i++) {
      int e = edges[i];
      const auto& xe = mesh_->edge_centroid(e);
      const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
      double length = mesh_->edge_length(e);

      for (auto jt = pe.begin(); jt < pe.end(); ++jt) {
        int m2 = jt.MonomialSetOrder();
        FunctionPower efunc(tau[k] / length, m2);

        funcs[0] = &cmono;
        funcs[1] = &efunc;

        N(row, col) = numi.IntegrateFunctionsEdge(e, funcs, m + m2) / length;
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
  std::vector<Basis_Regularized<SurfaceMiniMesh> > vbasisf;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

  Tensor Idf(d_ - 1, 2);
  Idf.MakeDiagonal(1.0);

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
    auto L2f = NN * NfT;
    vL2f.push_back(L2f);

    Basis_Regularized<SurfaceMiniMesh> basis_f;
    basis_f.Init(surf_mesh, f, order_ + 1, integrals_.poly());
    vbasisf.push_back(basis_f);

    NumericalIntegration<SurfaceMiniMesh> numi_f(surf_mesh);
    GrammMatrix(numi_f, order_ + 1, integrals_, basis_f, MG);
    vMGf.push_back(Idf ^ MG);
  }

  // -- L2-projector in cell
  Tensor Idc(d_, 2);
  Idc.MakeDiagonal(1.0);

  DenseMatrix MGc, sMGc, rMGc, NT, L2c;
  integrals_.set_id(c);
  GrammMatrix(numi, order_ + 1, integrals_, basis, MGc);
  sMGc = MGc.SubMatrix(0, ndc, 0, ndc);

  NT.Transpose(N);
  auto NN = NT * N;
  NN.InverseMoorePenrose();
  L2c = NN * NT;

  // -- curl matrix combined with L2 projector
  {
    int mdc = PolynomialSpaceDimension(d_, order_ - 1);
    auto C = Curl3DMatrix(d_, order_);
    auto Mtmp = MGc.SubMatrix(0, MGc.NumRows(), 0, mdc);
    rMGc = (Idc ^ Mtmp) * C * L2c;
  }

  // -----------------
  // assemble matrix R
  // ------------------
  DenseMatrix R(ndof_S, ndc * d_);
  R.PutScalar(0.0);

  for (auto it = it0; it < it1; ++it) {
    int k = it.VectorComponent();
    int m = it.PolynomialPosition();
    int col = it.VectorPolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Monomial q(d_, index, factor);
    q.set_origin(xc);

    // vector decomposition of vector (0, q, 0) with q at k-th position
    VectorPolynomial p1;
    Polynomial p2;
    VectorDecomposition3DCurl(q, k, p1, p2);

    // contributions from faces: int_f (p1^x^n . Pi_f(v))
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
      int nrows = vL2f[n].NumRows();
      int ncols = vL2f[n].NumCols();
      int krows = vMGf[n].NumRows();
      WhetStone::DenseVector w(krows), p0v(ncols);

      Polynomial xyzn(xyz[0]);
      for (int i = 0; i < d_; ++i) xyzn(i + 1) = normal[i];  // xyz * normal
      VectorPolynomial p3D = p1 * xyzn - xyz * (p1 * normal);
      VectorPolynomial p2D = ProjectVectorPolynomialOnManifold(p3D, xf, *vcoordsys[n]->tau());
      auto v = ExpandCoefficients(p2D);
      v /= area;

      int stride1 = v.NumRows() / (d_ - 1);
      int stride2 = PolynomialSpaceDimension(d_ - 1, order_ + 1);
      v.Regroup(stride1, stride2);

      // calculate one factor in the L2 inner product
      vbasisf[n].ChangeBasisNaturalToMy(v, d_ - 1);

      vMGf[n].Multiply(v, w, false);

      // reduce polynomial degree by one
      stride1 = nrows / (d_ - 1);
      w.Regroup(stride2, stride1);

      vL2f[n].Multiply(w, p0v, true);

      for (int i = 0; i < p0v.NumRows(); ++i) {
        R(map[i], col) += p0v(i) * fdirs[n];
      }
    }

    // first contribution from cell: int_c ((p2 x) . Pi_c(v))
    if (order_ > 0) {
      int nrows = L2c.NumRows(); 
      int ncols = L2c.NumCols(); 
      WhetStone::DenseVector w(nrows), p2v(ncols);

      VectorPolynomial p3D = xyz * p2;
      auto v = ExpandCoefficients(p3D);

      int stride1 = v.NumRows() / d_;
      int stride2 = ndc;
      v.Regroup(stride1, stride2);

      sMGc.BlockMultiply(v, w, false);
      L2c.Multiply(w, p2v, true);

      double factor = basis.monomial_scales()[1];
      for (int i = 0; i < ncols; ++i) {
        R(i, col) += p2v(i) / factor;
      }
    }

    // second contribution from cell: int_c (p1^x . curl Pi_c(v))
    if (order_ > 0) {
      int nrows = rMGc.NumRows(); 
      int ncols = rMGc.NumCols(); 
      WhetStone::DenseVector p1v(ncols);

      VectorPolynomial p3D = p1 ^ xyz;
      auto v = ExpandCoefficients(p3D);

      int stride1 = v.NumRows() / d_;
      int stride2 = nrows / d_;
      v.Regroup(stride1, stride2);

      rMGc.Multiply(v, p1v, true);

      double factor = basis.monomial_scales()[1];
      for (int i = 0; i < ncols; ++i) {
        R(i, col) += p1v(i) / factor;
      }
    }
  }

  // DenseMatrix X;
  // X.Transpose(N);
  // PrintMatrix(X * R, "%12.5f", X.NumRows());

  // calculate Mc = R (R^T N)^{-1} R^T 
  DenseMatrix RT;
  RT.Transpose(R);

  Tensor Kinv(K);
  Kinv.Inverse();

  sMGc.InverseSPD();
  Mc = R * ((Idc * Kinv) ^ sMGc) * RT;

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


/* ******************************************************************
* Projector on a mesh face.
****************************************************************** */
void VEM_NedelecSerendipityType2::ProjectorFace_(
    int f, const std::vector<VectorPolynomial>& ve,
    const ProjectorType type, const Polynomial* moments, VectorPolynomial& uf)
{
  const auto& xf = mesh_->face_centroid(f);
  const auto& normal = mesh_->face_normal(f);
  auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);

  Teuchos::RCP<const SurfaceMiniMesh> surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));

  std::vector<VectorPolynomial> vve;
  for (int i = 0; i < ve.size(); ++i) {
    auto tmp = ProjectVectorPolynomialOnManifold(ve[i], xf, *coordsys->tau());
    vve.push_back(tmp);
  }

  ProjectorCell_<SurfaceMiniMesh>(surf_mesh, f, vve, vve, type, moments, uf);
  uf.ChangeOrigin(AmanziGeometry::Point(d_ - 1));
  for (int i = 0; i < uf.NumRows(); ++i) {
    uf[i].InverseChangeCoordinates(xf, *coordsys->tau());  
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

