/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#include <cmath>
#include <tuple>
#include <vector>

// Amanzi
#include "MeshLight.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "Basis_Regularized.hh"
#include "GrammMatrix.hh"
#include "MFD3D_CrouzeixRaviartAnyOrder.hh"
#include "NumericalIntegration.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
MFD3D_CrouzeixRaviartAnyOrder::MFD3D_CrouzeixRaviartAnyOrder(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
  : MFD3D(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> MFD3D_CrouzeixRaviartAnyOrder::schema() const
{
  int nk = PolynomialSpaceDimension(d_ - 1, order_ - 1);
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::MOMENT, nk));

  if (order_ > 1) {
    nk = PolynomialSpaceDimension(d_, order_ - 2);
    items.push_back(std::make_tuple(AmanziMesh::CELL, DOF_Type::MOMENT, nk));
  }

  return items;
}


/* ******************************************************************
* High-order consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_CrouzeixRaviartAnyOrder::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c); 
  double volume = mesh_->cell_volume(c); 

  // calculate degrees of freedom 
  Polynomial poly(d_, order_), pf(d_ - 1, order_ - 1), pc;
  if (order_ > 1) {
    pc.Reshape(d_, order_ - 2);
  }
  int nd = poly.size();
  int ndf = pf.size();
  int ndc = pc.size();

  int ndof = nfaces * ndf + ndc;
  N.Reshape(ndof, nd);
  R_.Reshape(ndof, nd);
  Ac.Reshape(ndof, ndof);
  G_.Reshape(nd, nd);

  // select regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_, ptmp);

  // pre-calculate integrals of natural monomials 
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_ - 2, integrals_);

  // populate matrices N and R
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) { 
    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d_, index, factor);
    cmono.set_origin(xc);  

    // N and R: degrees of freedom on faces 
    auto grad = Gradient(cmono);
     
    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(0);

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      AmanziGeometry::Point normal = mesh_->face_normal(f);

      // local coordinate system with origin at face centroid
      SurfaceCoordinateSystem coordsys(xf, normal);
      const auto& tau = *coordsys.tau();
      normal *= dirs[i];

      auto conormal = K * normal;
      Polynomial tmp = grad * conormal;
      tmp.ChangeCoordinates(xf, tau);

      for (int n = 0; n < tmp.size(); ++n) {
        R_(row + n, col) = tmp(n);
      }

      for (auto jt = pf.begin(); jt < pf.end(); ++jt) {
        const int* jndex = jt.multi_index();
        Polynomial fmono(d_ - 1, jndex, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        int n = jt.PolynomialPosition();
        N(row + n, col) = numi.IntegratePolynomialsFace(f, polys) / area;
      }
      row += ndf;
    }

    // N and R: degrees of freedom in cells
    if (cmono.order() > 1) {
      VectorPolynomial Kgrad = K * grad;
      Polynomial tmp = Divergence(Kgrad);

      for (auto jt = tmp.begin(); jt < tmp.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int n = jt.PolynomialPosition();

        R_(row + n, col) = -tmp(n) / basis.monomial_scales()[m] * volume;
      }
    }

    if (order_ > 1) {
      for (auto jt = pc.begin(); jt < pc.end(); ++jt) {
        int n = jt.PolynomialPosition();
        const int* jndex = jt.multi_index();

        int nm(0);
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        int np = MonomialSetPosition(d_, multi_index); 
        factor = basis.monomial_scales()[it.MonomialSetOrder()] *
                 basis.monomial_scales()[jt.MonomialSetOrder()];
        N(row + n, col) = integrals_.poly()(nm, np) * factor / volume; 
      }
    }
  }

  // Gramm matrix for gradients of polynomials 
  G_.Multiply(N, R_, true);

  // calculate R inv(G) R^T
  DenseMatrix RG(ndof, nd), Rtmp(nd, ndof);

  // to invert degenerate matrix, we add and subtruct positive diagonal entry 
  G_(0, 0) = 1.0;
  G_.Inverse();
  G_(0, 0) = 0.0;
  RG.Multiply(R_, G_, false);

  Rtmp.Transpose(R_);
  Ac.Multiply(RG, Rtmp, false);

  return 0;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_CrouzeixRaviartAnyOrder::StiffnessMatrix(
    int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return 0;
}


/* ******************************************************************
* L2 projector of gradient on the space of polynomials of order k-1.
****************************************************************** */
void MFD3D_CrouzeixRaviartAnyOrder::ProjectorGradientCell_(
    int c, const std::vector<VectorPolynomial>& vf,
    const ProjectorType type, 
    const std::shared_ptr<DenseVector>& moments, MatrixPolynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  StiffnessMatrix(c, T, A);  

  // number of degrees of freedom
  Polynomial poly(d_, order_ -1);

  // create zero vector polynomial
  int dim = vf[0].size();
  uc.Reshape(d_, dim, d_, order_ - 1, true);

  std::vector<const PolynomialBase*> polys(2);
  NumericalIntegration numi(mesh_);

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_, ptmp);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < d_; ++j) {
      // calculate right-hand side matrix R for the L2 projector
      int md = poly.size();
      DenseVector v4(md), v5(md);

      v4.PutScalar(0.0);
      for (auto it = poly.begin(); it < poly.end(); ++it) {
        int row = it.PolynomialPosition();
        const int* index = it.multi_index();

        double factor = basis.monomial_scales()[it.MonomialSetOrder()];
        Polynomial cmono(d_, index, factor);
        cmono.set_origin(xc);  

        polys[0] = &cmono;

        // -- face contribution
        for (int n = 0; n < nfaces; ++n) {
          int f = faces[n];
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);
          double area = mesh_->face_area(f);

          polys[1] = &(vf[n][i]);
          double tmp = numi.IntegratePolynomialsFace(f, polys) / area;
          v4(row) += tmp * normal[j] * dirs[n];
        }

        // -- cell contribution
        if (order_ > 1) {
          auto grad = Gradient(cmono);

          for (auto jt = grad[j].begin(); jt < grad[j].end(); ++jt) {
            int m = jt.MonomialSetOrder();
            int s = jt.PolynomialPosition();
            v4(row) -= grad[j](s) / basis.monomial_scales()[m] * (*moments)(s) * volume;
          }
        }
      }

      // calculate coefficients of polynomial
      DenseMatrix M;
      GrammMatrix(poly, integrals_, basis, M);

      M.Inverse();
      M.Multiply(v4, v5, false);

      uc(i, j) = basis.CalculatePolynomial(mesh_, c, order_ - 1, v5);
    }
  }
}


/* ******************************************************************
* Degrees of freedom on face f.
****************************************************************** */
void MFD3D_CrouzeixRaviartAnyOrder::CalculateFaceDOFs_(
    int f, const Polynomial& vf, const Polynomial& pf,
    DenseVector& vdof, int& row)
{
  std::vector<const PolynomialBase*> polys(2);

  NumericalIntegration numi(mesh_);

  double area = mesh_->face_area(f);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);

  // local coordinate system with origin at face centroid
  SurfaceCoordinateSystem coordsys(xf, normal);

  polys[0] = &vf;

  for (auto it = pf.begin(); it < pf.end(); ++it) {
    const int* index = it.multi_index();
    Polynomial fmono(d_ - 1, index, 1.0);
    fmono.InverseChangeCoordinates(xf, *coordsys.tau());  

    polys[1] = &fmono;

    vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
    row++;
  }
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
void MFD3D_CrouzeixRaviartAnyOrder::ProjectorCell_(
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, Polynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  const auto& faces = mymesh->cell_get_faces(c);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mymesh->cell_centroid(c);
  double volume = mymesh->cell_volume(c);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  StiffnessMatrix(c, T, A);  

  // number of degrees of freedom
  Polynomial pf(d_ - 1, order_ - 1);
  int nd = G_.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  // calculate DOFs for boundary polynomial
  DenseVector vdof(ndof);

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mymesh, c, order_, ptmp);

  // populate matrices N and R
  int row(0);
  // degrees of freedom on faces
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    CalculateFaceDOFs_(f, vf[n], pf, vdof, row);
  }
  // degrees of freedom in cell
  if (ndof_c > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    AMANZI_ASSERT(ndof_c == v3.NumRows());

    for (int n = 0; n < ndof_c; ++n) {
      vdof(row + n) = v3(n);
    }
  }

  // calculate polynomial coefficients (in natural basis)
  DenseVector v4(nd), v5(nd);
  R_.Multiply(vdof, v4, true);
  G_.Multiply(v4, v5, false);

  uc = basis.CalculatePolynomial(mymesh, c, order_, v5);

  // uniqueness requires to specify constant in polynomial
  if (order_ == 1) {
    AmanziGeometry::Point grad(d_);
    for (int j = 0; j < d_; ++j) {
      grad[j] = uc(j + 1);
    }
    
    double a1(0.0), a2(0.0), tmp;
    for (int n = 0; n < nfaces; ++n) {  
      int f = faces[n];
      const AmanziGeometry::Point& xf = mymesh->face_centroid(f);
      double area = mymesh->face_area(f);
       
      tmp = vf[n].Value(xf) - grad * (xf - xc);
      a1 += tmp * area;
      a2 += area;
    }

    uc(0) = a1 / a2;
  } else if (order_ >= 2) {
    v4 = integrals_.poly().coefs();
    basis.ChangeBasisMyToNatural(v4);
    v4.Reshape(nd);
    uc(0) = vdof(row) - (v4 * v5) / volume;
  }

  // calculate L2 projector
  if (type == ProjectorType::L2 && ndof_c > 0) {
    v5(0) = uc(0);

    DenseMatrix M, M2;
    DenseVector v6(nd - ndof_c);
    Polynomial poly(d_, order_);
    NumericalIntegration numi(mymesh);

    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);
    GrammMatrix(poly, integrals_, basis, M);

    M2 = M.SubMatrix(ndof_c, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    const DenseVector& v3 = moments->coefs();
    for (int n = 0; n < ndof_c; ++n) {
      v4(n) = v3(n) * mymesh->cell_volume(c);
    }

    for (int n = 0; n < nd - ndof_c; ++n) {
      v4(ndof_c + n) = v6(n);
    }

    M.Inverse();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mymesh, c, order_, v5);
  }

  // set correct origin 
  uc.set_origin(xc);
}

}  // namespace WhetStone
}  // namespace Amanzi

