/*
  WhetStone, version 2.1
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
#include <vector>

// Amanzi
#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "Basis_Regularized.hh"
#include "CoordinateSystems.hh"
#include "DG_Modal.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_CrouzeixRaviart::H1consistencyLO_(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  N.Reshape(nfaces, d_ + 1);
  R_.Reshape(nfaces, d_ + 1);
  Ac.Reshape(nfaces, nfaces);

  // calculate matrix R
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    for (int k = 0; k < d_; k++) R_(n, k) = normal[k] * dirs[n];
    R_(n, d_) = 0.0;
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d_), v2(d_);
  double volume = mesh_->cell_volume(c);
  for (int n = 0; n < nfaces; ++n) {
    for (int k = 0; k < d_; k++) v1[k] = R_(n, k);
    v2 = K * v1;

    for (int m = n; m < nfaces; m++) {
      for (int k = 0; k < d_; k++) v1[k] = R_(m, k);
      Ac(n, m) = (v1 * v2) / volume;
    }
  }

  // calculate N
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    for (int k = 0; k < d_; k++) N(n, k) = xf[k] - xc[k];
    N(n, d_) = 1.0;  // additional column is added to the consistency condition
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix: the standard algorithm. 
****************************************************************** */
int MFD3D_CrouzeixRaviart::StiffnessMatrixLO_(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* High-order consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_CrouzeixRaviart::H1consistencyHO_(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
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
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_);

  // pre-calculate integrals of natural monomials 
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_ - 2, integrals_);

  // populate matrices N and R
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const Polynomial*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) { 
    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d_, index, factor);
    cmono.set_origin(xc);  

    // N and R: degrees of freedom on faces 
    VectorPolynomial grad;
    grad.Gradient(cmono);
     
    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(0);

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      AmanziGeometry::Point normal = mesh_->face_normal(f);
      FaceCoordinateSystem(normal, tau);
      normal *= dirs[i];

      auto conormal = K * normal;
      Polynomial tmp = grad * conormal;
      tmp.ChangeCoordinates(xf, tau);

      for (auto jt = tmp.begin(); jt < tmp.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int k = jt.MonomialSetPosition();
        int n = jt.PolynomialPosition();
        R_(row + n, col) = tmp(m, k);
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
        int k = jt.MonomialSetPosition();
        int n = jt.PolynomialPosition();

        R_(row + n, col) = -tmp(m, k) / basis.monomial_scales()[m] * volume;
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

        int np = poly.MonomialSetPosition(multi_index); 
        double factor = basis.monomial_scales()[it.MonomialSetOrder()] *
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

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_CrouzeixRaviart::StiffnessMatrixHO_(
    int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistencyHO_(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Energy projector on the space of linear polynomials in cell c.
****************************************************************** */
void MFD3D_CrouzeixRaviart::ProjectorCell_LO_(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double vol = mesh_->cell_volume(c);

  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, 1, true);
  }

  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    for (int i = 0; i < dim; ++i) {
      double tmp = vf[n][i].Value(xf) * dirs[n] / vol;

      for (int j = 0; j < d_; ++j) {
        uc[i](1, j) += tmp * normal[j];
      }
    }
  }

  // calculate projector's low-order term
  AmanziGeometry::Point grad(d_);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < d_; ++j) {
      grad[j] = uc[i](1, j);
    }
    
    double a1(0.0), a2(0.0), tmp;
    for (int n = 0; n < nfaces; ++n) {  
      int f = faces[n];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);
       
      tmp = vf[n][i].Value(xf) - grad * (xf - xc);
      a1 += tmp * area;
      a2 += area;
    }

    uc[i](0, 0) = a1 / a2;
  }

  // set the correct origin
  uc.set_origin(xc);
}


/* ******************************************************************
* Energy projector on space of linear polynomials in face f.
* Uniqueness requires to specify projector's value at face centroid.
****************************************************************** */
void MFD3D_CrouzeixRaviart::H1Face(
    int f, const AmanziGeometry::Point& p0,
    const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const
{
  Entity_ID_List edges;
  std::vector<int> dirs;

  mesh_->face_get_edges_and_dirs(f, &edges, &dirs);
  int nedges = edges.size();

  double area = mesh_->face_area(f);
  AmanziGeometry::Point fnormal = mesh_->face_normal(f);
  fnormal /= norm(fnormal);

  // create zero vector polynomial
  uf.resize(d_);
  for (int i = 0; i < d_; ++i) { 
    uf[i].Reshape(d_, 1, true);
  }

  AmanziGeometry::Point enormal(d_);

  for (int n = 0; n < nedges; ++n) {  
    int e = edges[n];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);

    enormal = tau^fnormal;

    for (int i = 0; i < d_; ++i) {
      double tmp = ve[n][i].Value(xe) * dirs[n] / area;

      for (int j = 0; j < d_; ++j) {
        uf[i](1, j) += tmp * enormal[j];
      }
    }
  }

  // fix the constant value
  const AmanziGeometry::Point& xf0 = mesh_->face_centroid(f);
  AmanziGeometry::Point zero(d_);

  for (int i = 0; i < d_; ++i) {
    uf[i](0, 0) = p0[i];
    uf[i].set_origin(xf0);
    uf[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
void MFD3D_CrouzeixRaviart::ProjectorCell_HO_(
    int c, const std::vector<VectorPolynomial>& vf,
    const Projectors::Type type,
    VectorPolynomial& moments, VectorPolynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  Entity_ID_List faces;
  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  StiffnessMatrixHO_(c, T, A);  

  // number of degrees of freedom
  Polynomial pf(d_ - 1, order_ - 1);
  int nd = G_.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, order_, true);
  }

  // calculate DOFs for boundary polynomial
  DenseVector vdof(ndof);
  NumericalIntegration numi(mesh_);

  // selecting regularized basis
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_);

  // populate matrices N and R
  for (int i = 0; i < dim; ++i) {
    int row(0);
    // degrees of freedom on faces
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      CalculateFaceDOFs_(f, vf[n][i], pf, vdof, row);
    }

    // degrees of freedom in cell
    if (ndof_c > 0) {
      DenseVector v3(ndof_c);
      moments[i].GetPolynomialCoefficients(v3);

      AMANZI_ASSERT(ndof_c == v3.NumRows());

      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = v3(n);
      }
    }

    // calculate polynomial coefficients (in natural basis)
    DenseVector v4(nd), v5(nd);
    R_.Multiply(vdof, v4, true);
    G_.Multiply(v4, v5, false);

    uc[i] = basis.CalculatePolynomial(mesh_, c, order_, v5);

    // uniqueness requires to specify constant in polynomial
    if (order_ == 1) {
      AmanziGeometry::Point grad(d_), zero(d_);
      for (int j = 0; j < d_; ++j) {
        grad[j] = uc[i](1, j);
      }
    
      double a1(0.0), a2(0.0), tmp;
      for (int n = 0; n < nfaces; ++n) {  
        int f = faces[n];
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        double area = mesh_->face_area(f);
       
        tmp = vf[n][i].Value(xf) - grad * (xf - xc);
        a1 += tmp * area;
        a2 += area;
      }

      uc[i](0, 0) = a1 / a2;
    } else if (order_ >= 2) {
      integrals_.poly().GetPolynomialCoefficients(v4);
      basis.ChangeBasisMyToNatural(v4);
      v4.Reshape(nd);
      uc[i](0, 0) = vdof(row) - (v4 * v5) / volume;
    }

    // calculate L2 projector
    if (type == Type::L2 && ndof_c > 0) {
      v5(0) = uc[i](0, 0);

      DG_Modal dg(order_, mesh_, "regularized");

      DenseMatrix M, M2;
      DenseVector v6(nd - ndof_c);
      dg.MassMatrix(c, T, integrals_, M);

      M2 = M.SubMatrix(ndof_c, nd, 0, nd);
      M2.Multiply(v5, v6, false);

      DenseVector v3(ndof_c);
      moments[i].GetPolynomialCoefficients(v3);

      for (int n = 0; n < ndof_c; ++n) {
        v4(n) = v3(n) * mesh_->cell_volume(c);
      }

      for (int n = 0; n < nd - ndof_c; ++n) {
        v4(ndof_c + n) = v6(n);
      }

      M.Inverse();
      M.Multiply(v4, v5, false);

      uc[i] = basis.CalculatePolynomial(mesh_, c, order_, v5);
    }

    // set correct origin 
    uc[i].set_origin(xc);
  }
}


/* ******************************************************************
* L2 projector of gradient on the space of polynomials of order k-1.
****************************************************************** */
void MFD3D_CrouzeixRaviart::ProjectorGradientCell_(
    int c, const std::vector<VectorPolynomial>& vf,
    const Projectors::Type type, 
    const std::shared_ptr<DenseVector>& moments, MatrixPolynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  StiffnessMatrixHO_(c, T, A);  

  // number of degrees of freedom
  Polynomial poly(d_, order_ -1), pf(d_ - 1, order_ - 1);
  int nd = G_.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].resize(d_);
    for (int j = 0; j < d_; ++j) { 
      uc[i][j].Reshape(d_, order_ - 1, true);
    }
  }

  std::vector<const Polynomial*> polys(2);
  NumericalIntegration numi(mesh_);

  // selecting regularized basis
  Basis_Regularized basis;
  basis.Init(mesh_, c, order_);

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
          VectorPolynomial grad(d_, d_);
          grad.Gradient(cmono);

          for (auto jt = grad[j].begin(); jt < grad[j].end(); ++jt) {
            int m = jt.MonomialSetOrder();
            int k = jt.MonomialSetPosition();
            int s = jt.PolynomialPosition();
            v4(row) -= grad[j](m, k) / basis.monomial_scales()[m] * (*moments)(s) * volume;
          }
        }
      }

      // calculate coefficients of polynomial
      DenseMatrix M;
      DG_Modal dg(order_ - 1, mesh_, "regularized");
      dg.MassMatrix(c, T, integrals_, M);

      M.Inverse();
      M.Multiply(v4, v5, false);

      uc[i][j] = basis.CalculatePolynomial(mesh_, c, order_ - 1, v5);
    }
  }
}


/* ******************************************************************
* Degrees of freedom on face f.
****************************************************************** */
void MFD3D_CrouzeixRaviart::CalculateFaceDOFs_(
    int f, const Polynomial& vf, const Polynomial& pf,
    DenseVector& vdof, int& row)
{
  std::vector<const Polynomial*> polys(2);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);

  NumericalIntegration numi(mesh_);

  const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
  double area = mesh_->face_area(f);

  // local coordinate system with origin at face centroid
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);
  FaceCoordinateSystem(normal, tau);

  polys[0] = &vf;

  for (auto it = pf.begin(); it < pf.end(); ++it) {
    const int* index = it.multi_index();
    Polynomial fmono(d_ - 1, index, 1.0);
    fmono.InverseChangeCoordinates(xf, tau);  

    polys[1] = &fmono;

    vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
    row++;
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

