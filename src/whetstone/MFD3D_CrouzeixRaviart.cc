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
#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "Basis_Regularized.hh"
#include "SurfaceCoordinateSystem.hh"
#include "GrammMatrix.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor parses the parameter list
****************************************************************** */
MFD3D_CrouzeixRaviart::MFD3D_CrouzeixRaviart(const Teuchos::ParameterList& plist,
                                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D(mesh),
    InnerProduct(mesh)
{
  order_ = plist.get<int>("method order");
}


/* ******************************************************************
* Consistency condition for stiffness matrix. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_CrouzeixRaviart::H1consistency(
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
int MFD3D_CrouzeixRaviart::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Energy projector on the space of linear polynomials in cell c.
****************************************************************** */
void MFD3D_CrouzeixRaviart::ProjectorCell_(
    int c, const std::vector<Polynomial>& vf, Polynomial& uc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double vol = mesh_->cell_volume(c);

  // create zero vector polynomial
  uc.Reshape(d_, 1, true);

  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    double tmp = vf[n].Value(xf) * dirs[n] / vol;

    for (int j = 0; j < d_; ++j) {
      uc(1, j) += tmp * normal[j];
    }
  }

  // calculate projector's low-order term
  AmanziGeometry::Point grad(d_);
  for (int j = 0; j < d_; ++j) {
    grad[j] = uc(1, j);
  }
    
  double a1(0.0), a2(0.0), tmp;
  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);
       
    tmp = vf[n].Value(xf) - grad * (xf - xc);
    a1 += tmp * area;
    a2 += area;
  }

  uc(0) = a1 / a2;

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

  for (int i = 0; i < d_; ++i) {
    uf[i](0) = p0[i];
    uf[i].set_origin(xf0);
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

