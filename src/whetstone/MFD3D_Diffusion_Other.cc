/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscaleneous discretization methods for diffusion.
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"
#include "WhetStoneMeshUtils.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* The conventional FV scheme for a general mesh.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseTPFA(int c, const Tensor& K, DenseMatrix& W)
{
  const auto& faces = mesh_->getCellFaces(c);
  const auto& dirs = mesh_->getCellFaceDirections(c);
  int nfaces = faces.size();

  W.Reshape(nfaces, nfaces);
  W.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  AmanziGeometry::Point a(d_);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f, c);
    //const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);

    a = xf - xc;
    // double s = mesh_->getFaceArea(f) * dirs[n] / norm(a);
    double s = mesh_->getFaceArea(f) / norm(a);
    double Knn = ((K * a) * normal) * s;
    double dxn = a * normal;
    W(n, n) = Knn / fabs(dxn);
  }

  return 0;
}


/* ******************************************************************
* The one-sided transmissibility coefficient. Any change to this 
* routine must be consistent with the above routine.
****************************************************************** */
double MFD3D_Diffusion::Transmissibility(int f, int c, const Tensor& K)
{
  int dir;
  const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
  AmanziGeometry::Point a(d_);

  const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
  const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f,  c, &dir);

  a = xf - xc;
  double s = mesh_->getFaceArea(f) * dir / norm(a);
  double Knn = ((K * a) * normal) * s;
  double dxn = a * normal;
  double W = Knn / fabs(dxn);

  return W;
}


/* ******************************************************************
* The debug version of the above FV scheme for a scalar tensor and
* an orthogonal brick element.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseDiagonal(int c, const Tensor& K, DenseMatrix& W)
{
  double volume = mesh_->getCellVolume(c);

  const auto& faces = mesh_->getCellFaces(c);
  int nfaces = faces.size();

  W.Reshape(nfaces, nfaces);
  W.PutScalar(0.0);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    double area = mesh_->getFaceArea(f);
    W(n, n) = nfaces * K(0, 0) * area * area / (d_ * volume);
  }
  return 0;
}


/* ******************************************************************
* Second-generation MFD method as inlemented in RC1.
****************************************************************** */
int MFD3D_Diffusion::MassMatrixInverseSO(int c, const Tensor& K, DenseMatrix& W)
{
  const auto& faces = mesh_->getCellFaces(c);
  const auto& fdirs = mesh_->getCellFaceDirections(c);
  int num_faces = faces.size();

  Entity_ID_List nodes, corner_faces;
  mesh_->getCellNodes(c, nodes);
  int nnodes = nodes.size();

  Tensor Kinv(K);
  Kinv.Inverse();

  // collect all corner matrices
  std::vector<Tensor> Mv;
  std::vector<double> cwgt;

  Tensor N(d_, 2), NK(d_, 2);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    node_get_cell_faces(*mesh_, v, c, Parallel_type::ALL, &corner_faces);
    int nfaces = corner_faces.size();
    if (nfaces < d_) {
      Errors::Message msg;
      msg << "WhetStone MFD3D_Diffusion: number of faces forming a corner is small.";
      Exceptions::amanzi_throw(msg);
    }

    for (int i = 0; i < d_; i++) {
      int f = corner_faces[i];
      N.SetColumn(i, mesh_->getFaceNormal(f));
    }
    double cwgt_tmp = fabs(N.Det());

    N.Inverse();
    NK = N * Kinv;

    N.Transpose();
    auto Mv_tmp = NK * N;
    Mv.push_back(Mv_tmp);

    for (int i = 0; i < d_; i++) {
      int f = corner_faces[i];
      cwgt_tmp /= mesh_->getFaceArea(f);
    }
    cwgt.push_back(cwgt_tmp);
  }

  // rescale corner weights
  double factor = 0.0;
  for (int n = 0; n < nnodes; n++) factor += cwgt[n];
  factor = mesh_->getCellVolume(c) / factor;

  for (int n = 0; n < nnodes; n++) cwgt[n] *= factor;

  // assemble corner matrices
  W.Reshape(num_faces, num_faces);
  W.PutScalar(0.0);
  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    node_get_cell_faces(*mesh_, v, c, Parallel_type::ALL, &corner_faces);

    Tensor& Mv_tmp = Mv[n];
    for (int i = 0; i < d_; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), corner_faces[i]));
      for (int j = i; j < d_; j++) {
        int l = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), corner_faces[j]));
        W(k, l) += Mv_tmp(i, j) * cwgt[n] * fdirs[k] * fdirs[l];
        W(l, k) = W(k, l);
      }
    }
  }
 
  // invert matrix W
  int ierr = W.Inverse();
  if (ierr != 0) {
    Errors::Message msg;
    msg << "WhetStone MFD3D_Diffusion: support operator generated bad elemental mass matrix.";
    Exceptions::amanzi_throw(msg);
  }

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi



