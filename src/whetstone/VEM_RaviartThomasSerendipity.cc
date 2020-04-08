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
#include "FunctionPower.hh"
#include "GrammMatrix.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "SurfaceMiniMesh.hh"
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
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D(mesh)
{
  order_ = plist.get<int>("method order");
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

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  double volume = mesh_->cell_volume(c);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  // selecting regularized basis (parameter integrals_ is not used)
  Basis_Regularized<AmanziMesh::Mesh> basis;
  basis.Init(mesh_, c, order_, integrals_.poly());

  NumericalIntegration<AmanziMesh::Mesh> numi(mesh_);
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
  std::vector<Basis_Regularized<SurfaceMiniMesh> > vbasisf;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

  MGf_.clear();
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
    auto surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));
    vcoordsys.push_back(coordsys);

    PolynomialOnMesh integrals_f;
    integrals_f.set_id(f);

    Basis_Regularized<SurfaceMiniMesh> basis_f;
    basis_f.Init(surf_mesh, f, order_ + 1, integrals_f.poly());
    vbasisf.push_back(basis_f);

    NumericalIntegration<SurfaceMiniMesh> numi_f(surf_mesh);
    GrammMatrix(numi_f, order_ + 1, integrals_f, basis_f, MG);
    vMGf.push_back(MG);

    auto S = MG.SubMatrix(0, ndf, 0, ndf);
    if (save_face_matrices_) MGf_.push_back(S);
    S.InverseSPD();
    vMGs.push_back(S);
  }

  // rows of matrix N are mean moments of gradients times normal
  // shift of a basis polynomial by a constant is not needed here 
  N.Reshape(ndof_S, d_ * ndc);

  // allocate sufficient memory to reduce memory copy inside vector reshape 
  int nrows = vMGf[0].NumRows();
  DenseVector v(nrows, 0.0), w(nrows, 0.0), u(ndf, 0.0);

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

      AmanziGeometry::Point knormal = K * normal;
      Polynomial poly(d_, index, factor * knormal[k]);
      poly.set_origin(xc);  

      poly.ChangeCoordinates(xf, *vcoordsys[i]->tau());  
      v = poly.coefs();
      vbasisf[i].ChangeBasisNaturalToMy(v);

      v.Reshape(nrows);
      vMGf[i].Multiply(v, w, false);

      // one scale area factor for normal, one for mean integral value
      double tmp = dirs[i] / area / area;
      for (auto l = 0; l < ndf; ++l) {
        N(row++, col) = tmp * w(l);
      }
    }
  }

  // L2 (serendipity-type)  projector 
  DenseMatrix NT;
  NT.Transpose(N);
  auto NN = NT * N;
  NN.InverseSPD();
  auto L2c = N * NN;

  // fixed vector
  DenseMatrix MGc;
  VectorPolynomial xyz(d_, d_, 1);
  if (order_ > 0) {
    for (int i = 0; i < d_; ++i) xyz[i](i + 1) = 1.0;
    xyz.set_origin(xc);

    integrals_.set_id(c);
    GrammMatrix(numi, order_, integrals_, basis, MGc);
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

  // calculate Mc = R (R^T N)^{-1} R^T 
  DenseMatrix RT;
  RT.Transpose(R);

  G_ = RT * N;
  G_.Inverse();

  Mc = R * G_ * RT;

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
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
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi

