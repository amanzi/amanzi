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
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "SurfaceMiniMesh.hh"
#include "Tensor.hh"
#include "VectorObjectsUtils.hh"
#include "VEM_NedelecSerendipityType2.hh"
#include "VEM_RaviartThomasSerendipity.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
VEM_NedelecSerendipityType2::VEM_NedelecSerendipityType2(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : DeRham_Edge(mesh),
    BilinearForm(mesh)
{
  // order of the maximum polynomial space
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> VEM_NedelecSerendipityType2::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::EDGE, DOF_Type::SCALAR, order_ + 1));

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
    return 0;
  }

  Entity_ID_List edges, faces, fedges;
  std::vector<int> fdirs, edirs;

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);

  // selecting regularized basis (parameter integrals_ is not used)
  Basis_Regularized<AmanziMesh::Mesh> basis;
  basis.Init(mesh_, c, order_ + 1, integrals_.poly());

  // calculate degrees of freedom: serendipity space S contains all boundary dofs
  Polynomial pc(d_, order_), pf(d_ - 1, order_);

  int ndc, nde;
  ndc = pc.size();
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
  MatrixOfDofs_(c, edges, basis, numi, N);

  // L2 projectors on faces
  std::vector<WhetStone::DenseMatrix> vL2f, vMGf;
  std::vector<Basis_Regularized<SurfaceMiniMesh> > vbasisf;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

  L2ProjectorsOnFaces_(c, K, faces, vL2f, vMGf, vbasisf, vcoordsys, order_ + 1);
  // MassMatricesOnFaces_(c, K, faces, vL2f, vMGf, vbasisf, vcoordsys, order_ + 1);

  // L2-projector in cell
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
    // contributions from faces: int_f (p1^x^n . v) 
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

      // method I
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

      // method II


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

      factor = basis.monomial_scales()[1];
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

      factor = basis.monomial_scales()[1];
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

  return 0;
}


/* ******************************************************************
* Mass matrix for edge-based discretization.
****************************************************************** */
void VEM_NedelecSerendipityType2::MatrixOfDofs_(
    int c, const Entity_ID_List& edges,
    const Basis_Regularized<AmanziMesh::Mesh>& basis,
    const NumericalIntegration<AmanziMesh::Mesh>& numi,
    DenseMatrix& N)
{
  int nedges = edges.size();
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  std::vector<double> moments;

  VectorPolynomialIterator it0(d_, d_, order_), it1(d_, d_, order_);
  it0.begin();
  it1.end();

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
      const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
      double length = mesh_->edge_length(e);

      numi.CalculatePolynomialMomentsEdge(e, cmono, order_, moments);
      for (int l = 0; l < moments.size(); ++l) {
        N(row, col) = moments[l] * tau[k] / length;
        row++;
      }
    }
  }
}


/* ******************************************************************
* Mass matrix for edge-based discretization.
****************************************************************** */
void VEM_NedelecSerendipityType2::L2ProjectorsOnFaces_(
    int c, const Tensor& K, const Entity_ID_List& faces,
    std::vector<WhetStone::DenseMatrix>& vL2f, 
    std::vector<WhetStone::DenseMatrix>& vMGf,
    std::vector<Basis_Regularized<SurfaceMiniMesh> >& vbasisf,
    std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> >& vcoordsys,
    int MGorder)
{
  int nfaces = faces.size();

  Tensor Idf(d_ - 1, 2);
  Idf.MakeDiagonal(1.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
    auto surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));
    vcoordsys.push_back(coordsys);

    DenseMatrix Nf, NfT, Mf, MG;
    L2consistency2D_<SurfaceMiniMesh>(surf_mesh, f, K, Nf, Mf, MG);
    if (MGorder == order_) vMGf.push_back(Idf ^ MG);

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
    if (MGorder == order_ + 1) vMGf.push_back(Idf ^ MG);
  }
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
  return 0;
}


/* ******************************************************************
* Stiffness matrix for edge-based discretization.
****************************************************************** */
int VEM_NedelecSerendipityType2::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix M, C;

  StiffnessMatrix(c, T, A, M, C);
  return 0;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. Curls in 2D and 3D are
* defined using exterior face normals.
****************************************************************** */
int VEM_NedelecSerendipityType2::StiffnessMatrix(
    int c, const Tensor& T, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C)
{
  Teuchos::ParameterList plist;
  plist.set<int>("method order", order_);

  VEM_RaviartThomasSerendipity rts(plist, mesh_);
  rts.MassMatrix(c, T, M);

  // populate curl matrix
  CurlMatrix(c, C);
  int ndofs_f = C.NumRows();
  int ndofs_e = C.NumCols();

  A.Reshape(ndofs_e, ndofs_e);
  DenseMatrix MC(ndofs_f, ndofs_e);

  MC.Multiply(M, C, false);
  A.Multiply(C, MC, true); 

  // rescaling (FIXME)
  Entity_ID_List faces, edges;
  mesh_->cell_get_faces(c, &faces);
  mesh_->cell_get_edges(c, &edges);

  int ndf = ndofs_f / faces.size();
  std::vector<double> areas;
  for (int i = 0; i < faces.size(); ++i) {
    double area = mesh_->face_area(faces[i]);
    for (int k = 0; k < ndf; ++k) areas.push_back(area);
  }

  for (int i = 0; i < ndofs_f; ++i) {
    for (int j = 0; j < ndofs_f; ++j) M(i, j) /= areas[i] * areas[j];
    for (int j = 0; j < ndofs_e; ++j) C(i, j) *= areas[i];
  }

  return 0;
}


/* ******************************************************************
* Curl matrix acts onto the space of fluxes.
****************************************************************** */
void VEM_NedelecSerendipityType2::CurlMatrix(int c, DenseMatrix& C)
{
  Entity_ID_List faces, nodes, fedges, edges;
  std::vector<int> fdirs, edirs, map;
  
  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  Polynomial pf(d_ - 1, order_);
  int ndf = pf.size();
  int nde = PolynomialSpaceDimension(d_ - 2, order_);
  int ncols = nedges * nde;
  int nrows = nfaces * ndf;

  C.Reshape(nrows, ncols);
  C.PutScalar(0.0);
 
  // precompute L2 projectors on faces
  std::vector<DenseMatrix> vL2f, vMGf;
  std::vector<Basis_Regularized<SurfaceMiniMesh> > vbasisf;
  std::vector<std::shared_ptr<SurfaceCoordinateSystem> > vcoordsys;

  Tensor K(d_, 1);
  K(0, 0) = 1.0;
  if (order_ > 0) {
    L2ProjectorsOnFaces_(c, K, faces, vL2f, vMGf, vbasisf, vcoordsys, order_);
  }

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);

    mesh_->face_to_cell_edge_map(f, c, &map);
    mesh_->face_get_edges_and_dirs(f, &fedges, &edirs);
    int nfedges = fedges.size();

    for (int m = 0; m < nfedges; ++m) {
      int e = fedges[m]; 
      double len = mesh_->edge_length(e);
      const auto& xe = mesh_->edge_centroid(e);

      int row = n * ndf;
      int col = map[m] * nde;

      // 0-th order moment
      C(row, col) = len * edirs[m] * fdirs[n] / area;

      // two 1-st order moments
      if (order_ > 0) {
        for (auto it = pf.begin(); it < pf.end(); ++it) {
          int pos = it.PolynomialPosition();
          if (pos == 0) continue;

          // surface terms
          double factor = vbasisf[n].monomial_scales()[it.MonomialSetOrder()];
          Polynomial fmono(d_ - 1, it.multi_index(), factor);
          std::vector<AmanziGeometry::Point> tau(1, vcoordsys[n]->Project(mesh_->edge_vector(e), false));
          fmono.ChangeCoordinates(vcoordsys[n]->Project(xe, true), tau);

          for (int k = 0; k < fmono.size(); ++k) {
            C(row + pos, col + k) += fmono(k) * len * edirs[m] * fdirs[n] / area;
          }
        }
      }
    }

    // volumetric term
    if (order_ > 0) {
      for (auto it = pf.begin(); it < pf.end(); ++it) {
        int pos = it.PolynomialPosition();
        double factor = vbasisf[n].monomial_scales()[it.MonomialSetOrder()];
        Polynomial fmono(d_ - 1, it.multi_index(), factor);

        auto rot = Rot2D(fmono);
        auto v1 = ExpandCoefficients(rot);

        int stride1 = v1.NumRows() / (d_ - 1);
        int stride2 = ndf;
        v1.Regroup(stride1, stride2);

        DenseVector v2(ndf * (d_ - 1)), v3(nfedges * nde);
        vMGf[n].Multiply(v1, v2, false);
        vL2f[n].Multiply(v2, v3, true);

        int l(0);
        for (int m = 0; m < nfedges; ++m) {
          int row = n * ndf;
          int col = map[m] * nde;

          for (int k = 0; k < nde; ++k) {
            C(row + pos, col + k) += v3(l) * fdirs[n] / area;
            l++;
          }
        }
      }
    }
  }
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

