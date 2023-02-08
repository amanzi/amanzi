/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Mimetic schemes for generalized polyhedra.
*/

#include "Mesh.hh"

#include "MFD3D_GeneralizedDiffusion.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem>
MFD3D_GeneralizedDiffusion::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::FACE, DOF_Type::SCALAR, d_));
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::CELL, DOF_Type::SCALAR, 1));
  return items;
}


/* ******************************************************************
* Consistency condition for inner product on a generized polyhedron.
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::L2consistency(int c,
                                          const Tensor& K,
                                          DenseMatrix& N,
                                          DenseMatrix& Mc,
                                          bool symmetry)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  int nx(d_ * nfaces);

  N.Reshape(nx, d_);
  Mc.Reshape(nx, nx);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  AmanziGeometry::Point v1(d_), v2(d_);
  AmanziMesh::Point_List vv(3), xm(3);

  // populate matrices R and N
  DenseMatrix R(N);
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double area = mesh_->getFaceArea(f);
    double area_div = norm(normal);

    CurvedFaceGeometry_(f, dirs[i], vv, xm);

    for (int k = 0; k < d_; ++k) {
      R(d_ * i, k) = area * xm[0][k] - area_div * xc[k];
      R(d_ * i + 1, k) = area * xm[1][k];
      R(d_ * i + 2, k) = area * xm[2][k];

      for (int l = 0; l < d_; ++l) { N(d_ * i + l, k) = vv[l][k]; }
    }
  }

  // upper triangular part of the consistency term
  Tensor Kinv(K);
  Kinv.Inverse();

  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < d_; ++k) v1[k] = R(i, k);
    v2 = Kinv * v1;

    for (int j = i; j < nx; ++j) {
      for (int k = 0; k < d_; ++k) v1[k] = R(j, k);
      Mc(i, j) = (v1 * v2) / volume;
    }
  }

  return 0;
}


/* ******************************************************************
* Mass matrix for genelized polyhedron
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistency(c, Kinv, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return 0;
}


/* ******************************************************************
* Mass matrix for genelized polyhedron
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::MassMatrixOptimized(int c, const Tensor& K, DenseMatrix& M)
{
  DenseMatrix N;

  Tensor Kinv(K);
  Kinv.Inverse();

  int ok = L2consistency(c, Kinv, N, M, true);
  if (ok) return ok;

  StabilityOptimized_(K, N, M);
  return 0;
}


/* ******************************************************************
* Consistency condition for inverse of inner product on a generized
* polyhedron.
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::L2consistencyInverse(int c,
                                                 const Tensor& K,
                                                 DenseMatrix& R,
                                                 DenseMatrix& Wc,
                                                 bool symmetry)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);

  int nfaces = faces.size();
  int nx(d_ * nfaces);

  R.Reshape(nx, d_);
  Wc.Reshape(nx, nx);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  double volume = mesh_->getCellVolume(c);

  AmanziGeometry::Point v1(d_), v2(d_);
  AmanziMesh::Point_List vv(3), xm(3);

  // populate matrices R and N
  DenseMatrix N(R);
  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    double area = mesh_->getFaceArea(f);
    double area_div = norm(normal);

    CurvedFaceGeometry_(f, dirs[i], vv, xm);

    for (int k = 0; k < d_; ++k) {
      R(d_ * i, k) = area * xm[0][k] - area_div * xc[k];
      R(d_ * i + 1, k) = area * xm[1][k];
      R(d_ * i + 2, k) = area * xm[2][k];

      for (int l = 0; l < d_; ++l) { N(d_ * i + l, k) = vv[l][k]; }
    }
  }

  // upper triangular part of the consistency term
  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = K * v1;

    for (int j = i; j < nx; ++j) {
      for (int k = 0; k < d_; ++k) v1[k] = N(j, k);
      Wc(i, j) = (v1 * v2) / volume;
    }
  }

  return 0;
}


/* ******************************************************************
* Inverse mass matrix for generalized polyhedron
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W)
{
  DenseMatrix R;

  int ok = L2consistencyInverse(c, K, R, W, true);
  if (ok) return ok;

  StabilityScalar_(R, W);
  return 0;
}


/* ******************************************************************
* Stiffness matrix is calculated by a hybridization algorithm.
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix M;
  MassMatrixInverse(c, K, M);

  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);

  int nfaces = faces.size();
  int nx(d_ * nfaces);

  // populate areas
  DenseVector area(nx), area_div(nx);
  area_div.PutScalar(0.0);

  for (int i = 0; i < nfaces; ++i) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    area_div(d_ * i) = norm(normal);

    double tmp = mesh_->getFaceArea(f);
    area(d_ * i) = tmp;
    area(d_ * i + 1) = tmp * dirs[i];
    area(d_ * i + 2) = tmp * dirs[i];
  }

  // populate stiffness matrix
  A.Reshape(nx + 1, nx + 1);

  double cntr(0.0);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nx; ++j) { A(i, j) = M(i, j) * area(i) * area(j); }

    double add(0.0);
    for (int j = 0; j < nx; ++j) { add -= M(i, j) * area_div(j); }
    A(nx, i) = A(i, nx) = add * area(i);

    cntr -= add * area_div(i);
  }
  A(nx, nx) = cntr;

  return 0;
}


/* ******************************************************************
* Divergence matrix.
****************************************************************** */
int
MFD3D_GeneralizedDiffusion::DivergenceMatrix(int c, DenseMatrix& A)
{
  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  A.Reshape(1, d_ * nfaces);
  A.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    A(0, d_ * n) = norm(normal) * dirs[n];
  }

  return 0;
}


/* ******************************************************************
* Geometry of a curved face
****************************************************************** */
void
MFD3D_GeneralizedDiffusion::CurvedFaceGeometry_(int f,
                                                int dirs,
                                                AmanziMesh::Point_List& vv,
                                                AmanziMesh::Point_List& xm)
{
  // local coordinate system uses external normal
  AmanziGeometry::Point normal(d_), xf(d_), v1(d_), v2(d_), v3(d_), p1(d_), p2(d_);

  normal = mesh_->getFaceNormal(f);
  normal /= norm(normal);

  if (fabs(normal[0]) > 0.1) {
    v1[0] = normal[1];
    v1[1] = -normal[0];
    v1[2] = 0.0;
  } else {
    v1[0] = 0.0;
    v1[1] = -normal[2];
    v1[2] = normal[1];
  }

  vv[0] = normal;
  vv[1] = v1 / norm(v1);
  vv[2] = normal ^ vv[1];

  vv[0] *= dirs; // exterior average normal

  // geometric center. We cannot use face_centroid
  auto nodes = mesh_->getFaceNodes(f);
  int nnodes = nodes.size();

  xf.set(0.0);
  for (int n = 0; n < nnodes; ++n) {
    p1 = mesh_->getNodeCoordinate(nodes[n]);
    xf += p1;
  }
  xf /= nnodes;

  // weighted center of gravity
  double area(0.0);

  xm[0].set(0.0, 0.0, 0.0);
  xm[1].set(0.0, 0.0, 0.0);
  xm[2].set(0.0, 0.0, 0.0);

  for (int n = 0; n < nnodes; ++n) {
    int m = (n + 1) % nnodes;

    p1 = mesh_->getNodeCoordinate(nodes[n]);
    p2 = mesh_->getNodeCoordinate(nodes[m]);

    v1 = xf - p1;
    v2 = xf - p2;
    v3 = v1 ^ v2;
    area += norm(v3) / 2;

    for (int k = 0; k < 3; ++k) {
      double s = dirs * (vv[k] * v3) / 2;
      xm[k] += s * (xf + p1 + p2) / 3;
    }
  }
  xm[0] /= area;
  xm[1] /= area;
  xm[2] /= area;
}

} // namespace WhetStone
} // namespace Amanzi
