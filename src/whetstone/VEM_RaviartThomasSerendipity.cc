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
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> VEM_RaviartThomasSerendipity::schema() const
{
  std::vector<SchemaItem> items;
  int nk = PolynomialSpaceDimension(d_ - 1, order_ - 1);
  items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, nk));

  return items;
}


/* ******************************************************************
* VEM scheme:: mass matrix for second-order scheme.
****************************************************************** */
int VEM_RaviartThomasSerendipity::L2consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
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
  numi.UpdateMonomialIntegralsCell(c, order_, integrals_);

  // calculate degrees of freedom: serendipity space S contains all boundary dofs
  Polynomial pc(d_, order_), pf(d_ - 1, order_ - 1);

  int ndc, ndf;
  ndc = pc.size();
  ndf = PolynomialSpaceDimension(d_ - 1, order_ - 1);
  int ndof_S = nfaces * ndf;

  // pre-compute Gramm matrices
  DenseMatrix MG;
  Polynomial qf(d_ - 1, order_);
  std::vector<WhetStone::DenseMatrix> vMGf, vMGs;
  std::vector<Basis_Regularized<SurfaceMiniMesh> > vbasisf;
  std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> > vcoordsys;

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
    basis_f.Init(surf_mesh, f, order_, integrals_f.poly());
    vbasisf.push_back(basis_f);

    NumericalIntegration<SurfaceMiniMesh> numi_f(surf_mesh);
    GrammMatrix(numi_f, order_, integrals_f, basis_f, MG);
    vMGf.push_back(MG);
 
    auto S = MG.SubMatrix(0, ndf, 0, ndf);
    S.InverseSPD();
    vMGs.push_back(S);
  }

  // rows of matrix N are mean moments of gradients times normal
  // shift of a basis polynomial by a constant is not needed here 
  N.Reshape(ndof_S, ndc - 1);

  // allocate sufficient memory to reduce memory copy inside vector reshape 
  int nrows = vMGf[0].NumRows();
  DenseVector v(nrows), w(nrows), u(ndf);

  PolynomialIterator it0 = pc.begin();
  ++it0;

  for (auto it = it0; it < pc.end(); ++it) { 
    int m = it.MonomialSetOrder();
    int col = it.PolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[m];
    Monomial cmono(d_, index, factor);
    cmono.set_origin(xc);  

    auto grad = Gradient(cmono);

    int row(0);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      AmanziGeometry::Point knormal = K * normal;
      auto poly = grad * knormal;
      poly.ChangeCoordinates(xf, *vcoordsys[i]->tau());  
      v = poly.coefs();
      vbasisf[i].ChangeBasisNaturalToMy(v);

      v.Reshape(nrows);
      vMGf[i].Multiply(v, w, false);

      // one scale area factor for normal, one for mean integral value
      double tmp = dirs[i] / area / area;
      for (auto k = 0; k < ndf; ++k) {
        N(row++, col - 1) = tmp * w(k);
      }
    }
  }

  // assemble matrix R
  // for 2nd-order scheme, basis polynomial is shifted by an orthogonalization
  // constant to annihilate the volume integral
  DenseMatrix R(ndof_S, ndc - 1);

  for (auto it = it0; it < pc.end(); ++it) {
    int m = it.MonomialSetOrder();
    int col = it.PolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[m];
    Monomial cmono(d_, index, factor);
    cmono.set_origin(xc);

    int row(0);
    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];
      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

      Polynomial poly(cmono);
      poly(0) -= (integrals_.poly())(col) / volume;
      poly.ChangeCoordinates(xf, *vcoordsys[i]->tau());  
      v = poly.coefs();
      vbasisf[i].ChangeBasisNaturalToMy(v);

      v.Reshape(nrows);
      w.Reshape(nrows);
      vMGf[i].Multiply(v, w, false);
      
      w.Reshape(ndf);
      vMGs[i].Multiply(w, u, false);

      for (int k = 0; k < ndf; ++k) {
        R(row++, col - 1) = u(k) * dirs[i] * area; 
      }
    }
  }
std::cout << R << " " << vMGf[0] << std::endl;

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

