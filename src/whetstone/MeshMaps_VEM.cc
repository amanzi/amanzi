/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Maps between mesh objects located on different meshes, e.g. two states 
  of a deformable mesh: virtual element implementation.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMaps_VEM.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity in cell c
****************************************************************** */
void MeshMaps_VEM::VelocityCell(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& vc) const
{
  vc.resize(d_);
  for (int i = 0; i < d_; ++i) vc[i].Reshape(d_, 1);
  
  Entity_ID_List nodes;
  mesh0_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  AmanziGeometry::Point px;
  std::vector<AmanziGeometry::Point> x1, x2, u;
  for (int i = 0; i < nnodes; ++i) {
    mesh0_->node_get_coordinates(nodes[i], &px);
    x1.push_back(px);
  }
  for (int i = 0; i < nnodes; ++i) {
    mesh1_->node_get_coordinates(nodes[i], &px);
    x2.push_back(px);
  }

  // calculate velocity u(X) = F(X) - X
  LeastSquareFit(1, x1, x2, vc);

  for (int i = 0; i < d_; ++i) {
    vc[i](1, i) -= 1.0;
  }

  // new method for velocity calculation
  // VectorPolynomial tmp; 
  EllipticProjectorP1(c, vf, vc);
  // std::cout << c << " " << mesh0_->cell_centroid(c) << "\n" << vc[0] << vc[1] << tmp[0] << tmp[1] << "\n\n";
}


/* ******************************************************************
* Transformation of normal is defined completely by face data.
* NOTE: Limited to P1 elements.
****************************************************************** */
void MeshMaps_VEM::NansonFormula(
    int f, double t, const VectorPolynomial& v, VectorPolynomial& cn) const
{
  AmanziGeometry::Point p(d_);
  WhetStone::Tensor J(d_, 2);

  JacobianFaceValue(f, v, p, J);
  J *= t;
  J += 1.0 - t;
  p = J * mesh0_->face_normal(f);

  cn.resize(d_);
  for (int i = 0; i < d_; ++i) {
    cn[i].Reshape(d_, 0);
    cn[i](0, 0) = p[i];
  }
}


/* ******************************************************************
* Calculation of matrix of cofactors
****************************************************************** */
void MeshMaps_VEM::Cofactors(
    int c, double t, const VectorPolynomial& vc, MatrixPolynomial& C) const
{
  C.resize(2);
  for (int i = 0; i < d_; ++i) {
    C[i].resize(2);
    for (int j = 0; j < d_; ++j) {
      double sgn = (i == j) ? t : -t;
      C[i][j].Reshape(d_, 0);
      C[i][j](0, 0) = sgn * vc[1 - i](1, 1 - j);
    }
    C[i][i](0, 0) += 1.0;
  }
}


/* ******************************************************************
* Calculation of Jacobian.
****************************************************************** */
void MeshMaps_VEM::JacobianCellValue(
    int c, double t, const AmanziGeometry::Point& xref, Tensor& J) const
{
}


/* ******************************************************************
* Calculate determinant of a Jacobian. A prototype for the future 
* projection scheme. Currently, we return a number.
****************************************************************** */
void MeshMaps_VEM::JacobianDet(
    int c, double t, const std::vector<VectorPolynomial>& vf, Polynomial& vc) const
{
  AmanziGeometry::Point x(d_), cn(d_);
  WhetStone::Tensor J(d_, 2); 

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh0_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  double sum(0.0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    // calculate j J^{-t} N dA
    JacobianFaceValue(f, vf[n], x, J);

    J *= t;
    J += 1.0 - t;

    Tensor C = J.Cofactors();
    cn = C * mesh0_->face_normal(f); 

    const AmanziGeometry::Point& xf0 = mesh0_->face_centroid(f);
    const AmanziGeometry::Point& xf1 = mesh1_->face_centroid(f);
    sum += (xf0 + t * (xf1 - xf0)) * cn * dirs[n];
  }
  sum /= 2 * mesh0_->cell_volume(c);

  vc.Reshape(d_, 0);
  vc(0, 0) = sum;
}


/* ******************************************************************
* Calculate Jacobian at point x of face f 
****************************************************************** */
void MeshMaps_VEM::JacobianFaceValue(
    int f, const VectorPolynomial& v,
    const AmanziGeometry::Point& x, Tensor& J) const
{
  // FIXME x is not used
  for (int i = 0; i < d_; ++i) {
    for (int j = 0; j < d_; ++j) {
      J(i, j) = v[i](1, j);
    }
  }
  J += 1.0;
}

}  // namespace WhetStone
}  // namespace Amanzi

