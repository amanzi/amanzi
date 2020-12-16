/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  High-order 3D serendipity Raviart Thomas element: degrees of freedom
  are moments on faces and inside cell. Degrees of freedom are ordered
  as follows:
    (1) moments on faces
    (2) moments inside cell
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
#include "SingleFaceMesh.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"
#include "VectorObjectsUtils.hh"
#include "VEM_RaviartThomasSerendipity.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
VEM_RaviartThomasSerendipity::VEM_RaviartThomasSerendipity(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
  : BilinearForm(mesh),
    save_face_matrices_(false)
{
  order_ = plist.get<int>("method order");
  if (plist.isParameter("save face matrices"))
    save_face_matrices_ = plist.get<bool>("save face matrices");
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> VEM_RaviartThomasSerendipity::schema() const
{
  std::vector<SchemaItem> items;
  int nk = PolynomialSpaceDimension(d_ - 1, order_);
  items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, nk));

  return items;
}


/* ******************************************************************
* VEM scheme:: mass matrix for second-order scheme.
****************************************************************** */
int VEM_RaviartThomasSerendipity::L2consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  // div v is no longer a constant for order > 1
  AMANZI_ASSERT(order_ < 2);

  Tensor Kinv(K);
  Kinv.Inverse();

  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  double volume = mesh_->cell_volume(c);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  // selecting regularized basis (parameter integrals_ is not used)
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_, integrals_.poly());

  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, order_ + 1, integrals_);

  // calculate degrees of freedom: serendipity space S contains all boundary dofs
  Polynomial pc(d_, order_), pf(d_ - 1, order_);

  int ndc, ndf;
  ndc = pc.size();
  ndf = PolynomialSpaceDimension(d_ - 1, order_);
  int ndof_S = nfaces * ndf;

  // pre-compute Gramm matrices
  DenseMatrix MG;
  std::vector<WhetStone::DenseMatrix> vMGf, vMGs;
  std::vector<Basis_Regularized> vbasisf;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

  MGf_.clear();
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
    auto surf_mesh = Teuchos::rcp(new SingleFaceMesh(mesh_, f, *coordsys));
    vcoordsys.push_back(coordsys);

    PolynomialOnMesh integrals_f;
    integrals_f.set_id(0);  // this is a one-cell mesh

    Basis_Regularized basis_f;
    basis_f.Init(surf_mesh, 0, order_ + 1, integrals_f.poly());
    vbasisf.push_back(basis_f);

    NumericalIntegration numi_f(surf_mesh);
    GrammMatrix(0, numi_f, order_ + 1, integrals_f, basis_f, MG);
    vMGf.push_back(MG);

    auto S = MG.SubMatrix(0, ndf, 0, ndf);
    if (save_face_matrices_) MGf_.push_back(S);
    S.InverseSPD();
    vMGs.push_back(S);
  }

  // allocate sufficient memory to reduce memory copy inside vector reshape 
  int nrows = vMGf[0].NumRows();
  DenseVector v(nrows, 0.0), w(nrows, 0.0), u(ndf, 0.0);

  // iterators
  VectorPolynomialIterator it0(d_, d_, order_), it1(d_, d_, order_);
  it0.begin();
  it1.end();

  // rows of matrix N are mean moments of gradients times normal
  // shift of a basis polynomial by a constant is not needed here 
  N.Reshape(ndof_S, d_ * ndc);
  ComputeN_(c, faces, dirs, basis, vbasisf, vcoordsys, vMGf, Kinv, N);

  Tensor Idc(d_, 1);
  Idc(0, 0) = 1.0;

  DenseMatrix N1(ndof_S, d_ * ndc);
  ComputeN_(c, faces, dirs, basis, vbasisf, vcoordsys, vMGf, Idc, N1);

  // L2 (serendipity-type)  projector 
  DenseMatrix N1T;
  N1T.Transpose(N1);
  auto NN = N1T * N1;
  NN.InverseSPD();
  auto L2c = N1 * NN;

  // fixed vector
  DenseMatrix MGc;
  VectorPolynomial xyz(d_, d_, 1);
  if (order_ > 0) {
    for (int i = 0; i < d_; ++i) xyz[i](i + 1) = 1.0;
    xyz.set_origin(xc);

    GrammMatrix(c, numi, order_, integrals_, basis, MGc);
  }

  // assemble matrix R
  Polynomial p1;
  VectorPolynomial p2;

  DenseMatrix R(ndof_S, d_ * ndc);

  for (auto it = it0; it < it1; ++it) {
    int k = it.VectorComponent();
    int m = it.MonomialSetOrder();
    int col = it.VectorPolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[m];
    Monomial cmono(d_, index, factor);
    cmono.set_origin(xc);

    // decompose vector monomial and othogonalize p1 (a monomial)
    int pos = VectorDecomposition3DGrad(cmono, k, p1, p2);
    p1(0) -= p1(pos) * (integrals_.poly())(pos) / volume;

    int row(0);
    // contribution from faces: int_f {p2 (Pi_f(v) . n)}
    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];
      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

      Polynomial poly(p1);
      poly.ChangeCoordinates(xf, *vcoordsys[i]->tau());  
      v = poly.coefs();
      vbasisf[i].ChangeBasisNaturalToMy(v);

      v.Reshape(nrows);
      w.Reshape(nrows);
      vMGf[i].Multiply(v, w, false);
      
      w.Reshape(ndf);
      vMGs[i].Multiply(w, u, false);

      for (int l = 0; l < ndf; ++l) {
        R(row++, col) = u(l) * area;
      }
    }

    // contribution from cell: int_c {(x ^ p2) . Pi_c(v)}
    if (m > 0) {
      WhetStone::DenseVector w3(d_ * ndc), v3(ndof_S);

      VectorPolynomial p3D = xyz ^ p2;
      auto u3 = ExpandCoefficients(p3D);

      int stride1 = u3.NumRows() / d_;
      int stride2 = ndc;
      basis.ChangeBasisNaturalToMy(u3, d_);
      u3.Regroup(stride1, stride2);

      MGc.BlockMultiply(u3, w3, false);
      L2c.Multiply(w3, v3, false);

      for (int i = 0; i < ndof_S; ++i) {
        R(i, col) += v3(i);
      }
    }
  }

  // DenseMatrix X;
  // X.Transpose(N);
  // PrintMatrix(X * R, "%12.8f", X.NumRows());
  // GrammMatrix(numi, order_, integrals_, basis, G_);
  // PrintMatrix(G_, "%12.8f", G_.NumRows());

  DenseMatrix RT;
  RT.Transpose(R);

  G_ = RT * N;
  G_.Inverse();

  Mc = R * G_ * RT;

  return 0;
}


/* ******************************************************************
* Mass matrix for edge-based discretization.
****************************************************************** */
int VEM_RaviartThomasSerendipity::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, K, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
void VEM_RaviartThomasSerendipity::ProjectorCell_(
    int c, const std::vector<VectorPolynomial>& ve,
    const std::vector<VectorPolynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, VectorPolynomial& uc)
{
  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_, ptmp);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  T(0, 0) = 1.0;

  DenseMatrix N, Mc;
  L2consistency(c, T, N, Mc, true);

  // degrees of freedom: serendipity space S contains all boundary dofs
  // plus a few internal dofs that depend on the value of eta.
  int ndof = N.NumRows();
  int ndof_cs = 0;  // required cell moments
  int ndof_s(ndof);  // serendipity dofs

  // extract submatrix
  int ncols = N.NumCols();
  DenseMatrix Ns, NN(ncols, ncols);
  Ns = N.SubMatrix(0, ndof_s, 0, ncols);

  NN.Multiply(Ns, Ns, true);
  NN.InverseSPD();

  // calculate degrees of freedom (Ns^T Ns)^{-1} Ns^T v
  // for consistency with other code, we use v5 for polynomial coefficients
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  DenseVector v1(ncols), v5(ncols);

  DenseVector vdof(ndof_s + ndof_cs);
  CalculateDOFsOnBoundary(c, ve, vf, vdof);

  // DOFs inside cell: copy moments from input data
  if (ndof_cs > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    AMANZI_ASSERT(ndof_cs == v3.NumRows());

    for (int n = 0; n < ndof_cs; ++n) {
      vdof(ndof_s + n) = v3(n);
    }
  }

  Ns.Multiply(vdof, v1, true);
  NN.Multiply(v1, v5, false);

  // this gives the least square projector
  int stride = v5.NumRows() / d_;
  DenseVector v4(stride);

  uc.resize(d_);
  for (int k = 0; k < d_; ++k) {
    for (int i = 0; i < stride; ++i) v4(i) = v5(k * stride + i);
    uc[k] = basis.CalculatePolynomial(mesh_, c, order_, v4);
  }

  // set correct origin 
  uc.set_origin(xc);
}


/* ******************************************************************
* Calculate boundary degrees of freedom in 2D and 3D.
****************************************************************** */
void VEM_RaviartThomasSerendipity::CalculateDOFsOnBoundary(
    int c, const std::vector<VectorPolynomial>& ve,
    const std::vector<VectorPolynomial>& vf, DenseVector& vdof)
{
  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  std::vector<const WhetStoneFunction*> funcs(2);
  NumericalIntegration numi(mesh_);

  // number of moments on faces
  std::vector<double> moments;

  int row(0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);
    const auto& normal = mesh_->face_normal(f);

    auto poly = vf[n] * normal;

    numi.CalculatePolynomialMomentsFace(f, poly, order_, moments);
    for (int k = 0; k < moments.size(); ++k) {
      vdof(row) = dirs[n] * moments[k] / area;
      row++;
    }
  }
}


/* ******************************************************************
* Mass matrix for face-based discretization.
****************************************************************** */
void VEM_RaviartThomasSerendipity::ComputeN_(
    int c, const Entity_ID_List& faces, const std::vector<int>& dirs,
    const Basis_Regularized& basis,
    const std::vector<Basis_Regularized>& vbasisf,
    const std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> >& vcoordsys,
    const std::vector<WhetStone::DenseMatrix>& vMGf,
    const Tensor& Kinv, DenseMatrix& N)
{
  int nfaces = faces.size();
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  int ndf = PolynomialSpaceDimension(d_ - 1, order_);

  // allocate sufficient memory to reduce memory copy inside vector reshape 
  int nrows = vMGf[0].NumRows();
  DenseVector v(nrows, 0.0), w(nrows, 0.0);

  // iterators
  VectorPolynomialIterator it0(d_, d_, order_), it1(d_, d_, order_);
  it0.begin();
  it1.end();

  for (auto it = it0; it < it1; ++it) { 
    int k = it.VectorComponent();
    int m = it.MonomialSetOrder();
    int col = it.VectorPolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[m];

    int row(0);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      AmanziGeometry::Point knormal = Kinv * normal;
      Polynomial poly(d_, index, factor * knormal[k]);
      poly.set_origin(xc);  

      poly.ChangeCoordinates(xf, *vcoordsys[i]->tau());  
      v = poly.coefs();
      vbasisf[i].ChangeBasisNaturalToMy(v);

      v.Reshape(nrows);
      vMGf[i].Multiply(v, w, false);

      // one factor "area" for normal, one for mean integral value
      double tmp = dirs[i] / area / area;
      for (auto l = 0; l < ndf; ++l) {
        N(row++, col) = tmp * w(l);
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

