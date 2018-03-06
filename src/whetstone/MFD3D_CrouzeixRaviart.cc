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
#include "CoordinateSystems.hh"
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
  Entity_ID_List edges, faces, face_edges;
  std::vector<int> dirs, map;

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  N.Reshape(nedges, d_ + 1);
  Ac.Reshape(nedges, nedges);

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int num_faces = faces.size();

  // to calculate matrix R, we use temporary matrix N 
  N.PutScalar(0.0);

  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    if (d_ == 2) {
      for (int k = 0; k < d_; k++) N(i, k) += normal[k] * dirs[i];
    } else {
      mesh_->face_to_cell_edge_map(f, c, &map);
      int num_face_edges = map.size();

      for (int j = 0; j < num_face_edges; j++) {
        for (int k = 0; k < d_; k++) N(map[j], k) += normal[k] * dirs[i];
      }
    }
  }

  // calculate R K R^T / volume
  AmanziGeometry::Point v1(d_), v2(d_);
  double volume = mesh_->cell_volume(c);
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; k++) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  // calculate N
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    for (int k = 0; k < d_; k++) N(i, k) = xe[k] - xc[k];
    N(i, d_) = 1;  // additional column is added to the consistency condition
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

  // pre-calculate integrals of monomials 
  NumericalIntegration numi(mesh_);
  integrals_.Reshape(d_, 2 * order_ - 2, true);

  for (int k = 0; k <= 2 * order_ - 2; ++k) {
    numi.IntegrateMonomialsCell(c, integrals_.monomials(k));
  }

  // populate matrices N and R
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const Polynomial*> polys(2);

  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) { 
    const int* index = it.multi_index();
    double factor = numi.MonomialNaturalScale(it.MonomialOrder(), volume);
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

      Polynomial tmp = grad * normal;
      tmp.ChangeCoordinates(xf, tau);

      for (auto jt = tmp.begin(); jt.end() <= tmp.end(); ++jt) {
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        int n = jt.PolynomialPosition();
        R_(row + n, col) = tmp(m, k);
      }

      for (auto jt = pf.begin(); jt.end() <= pf.end(); ++jt) {
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
      Polynomial tmp = cmono.Laplacian();
      for (auto jt = tmp.begin(); jt.end() <= tmp.end(); ++jt) {
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        int n = jt.PolynomialPosition();

        R_(row + n, col) = -tmp(m, k) * volume;
      }
    }

    if (order_ > 1) {
      for (auto jt = pc.begin(); jt.end() <= pc.end(); ++jt) {
        int n = jt.PolynomialPosition();
        const int* jndex = jt.multi_index();

        int nm(0);
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        // FIXME: use naturally scaled monomials for internal DOF
        double s = numi.MonomialNaturalScale(jt.MonomialOrder(), volume);

        const auto& coefs = integrals_.monomials(nm).coefs();
        N(row + n, col) = coefs[poly.MonomialPosition(multi_index)] / (volume * s); 
      }
    }
  }

  // set the Gramm-Schidt matrix for gradients of polynomials
  G_.PutScalar(0.0);

  // -- gradient of a naturally scaled polynomial needs correction
  double scale = numi.MonomialNaturalScale(1, volume);
   
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt.end() <= poly.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      
      int n(0);
      int multi_index[3];
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = index[i] + jndex[i];
        n += multi_index[i];
      }

      double sum(0.0), tmp;
      for (int i = 0; i < d_; ++i) {
        if (index[i] > 0 && jndex[i] > 0) {
          multi_index[i] -= 2;
          const auto& coefs = integrals_.monomials(n - 2).coefs();
          tmp = coefs[poly.MonomialPosition(multi_index)]; 
          sum += tmp * index[i] * jndex[i];
          multi_index[i] += 2;
        }
      }

      G_(l, k) = G_(k, l) = K(0, 0) * sum * scale * scale; 
    }
  }

  // calculate R inv(G) R^T
  DenseMatrix RG(ndof, nd), Rtmp(nd, ndof);

  // to invert generate matrix, we add and subtruct positive number
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
        uc[i].monomials(1).coefs()[j] += tmp * normal[j];
      }
    }
  }

  // calculate projector's low-order term
  AmanziGeometry::Point grad(d_), zero(d_);
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

  // clean-up: fix the origin
  for (int i = 0; i < dim; ++i) {
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Energy projector on space of linear polynomials in face f.
* Uniqueness requires to specify projector's value at face centroid.
****************************************************************** */
void MFD3D_CrouzeixRaviart::H1FaceHarmonic(
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
* Energy projector on space of polynomials of order k in cell c.
* Uniqueness require to specify its value at cell centroid.
****************************************************************** */
void MFD3D_CrouzeixRaviart::ProjectorCell_HO_(
    int c, const std::vector<VectorPolynomial>& vf,
    const Projectors::Type type, bool is_harmonic, 
    const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc)
{
  ASSERT(d_ == 2);

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
  Polynomial pf(d_ - 1, order_ - 1);
  int nd = G_.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  DenseMatrix Acf, Acc;
  if (ndof_c > 0 && is_harmonic) {
    Acf = A.SubMatrix(ndof_f, ndof, 0, ndof_f);
    Acc = A.SubMatrix(ndof_f, ndof, ndof_f, ndof);
    Acc.Inverse();
  }
  
  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, order_, true);
  }

  // calculate DOFs for boundary polynomial
  DenseVector vdof(ndof);
  std::vector<const Polynomial*> polys(2);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  NumericalIntegration numi(mesh_);

  for (int i = 0; i < dim; ++i) {
    int row(0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      FaceCoordinateSystem(normal, tau);

      polys[0] = &(vf[n][i]);

      for (auto it = pf.begin(); it.end() <= pf.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }

    // harmonic extension inside cell
    if (ndof_c > 0 && is_harmonic) {
      DenseVector v1(ndof_f), v2(ndof_c), v3(ndof_c);

      for (int n = 0; n < ndof_f; ++n) {
        v1(n) = vdof(n);
      }

      Acf.Multiply(v1, v2, false);
      Acc.Multiply(v2, v3, false);

      moments->Reshape(ndof_c);
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = -v3(n);
        (*moments)(n) = -v3(n);
      }
    }
    else if (ndof_c > 0 && !is_harmonic) {
      ASSERT(ndof_c == moments->NumRows());
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = (*moments)(n);
      }
    }

    // calculate polynomial coefficients
    DenseVector v4(nd), v5(nd);
    R_.Multiply(vdof, v4, true);
    G_.Multiply(v4, v5, false);

    uc[i].SetPolynomialCoefficients(v5);
    numi.ChangeBasisNaturalToRegular(c, uc[i]);

    // calculate the constant value
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
      integrals_.GetPolynomialCoefficients(v4);
      v4.Reshape(nd);
      uc[i](0, 0) = vdof(row) - (v4 * v5) / volume;
    }

    // change origin from centroid to zero
    AmanziGeometry::Point zero(d_);
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

