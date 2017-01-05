/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d_diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* The stable discretization for Stokes: vectors at nodes plus normal
* components at faces.
****************************************************************** */
int MFD3D_Diffusion::H1consistencyNodeFace(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int d = mesh_->space_dimension();
  AmanziGeometry::Point xv(d), v1(d);

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  // to calculate matrix R, we use temporary matrix N 
  N.PutScalar(0.0);

  Tensor T0(d, 2), T1(d, 2), T2(d, 2);
  T0(0, 0) = 1.0;
  T1(0, 1) = T1(1, 0) = 1.0;
  T2(1, 1) = 1.0;

  Tensor KT0(K * T0);
  Tensor KT1(K * T1);
  Tensor KT2(K * T2);

  const AmanziGeometry::Point& xm = mesh_->cell_centroid(c);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    const AmanziGeometry::Point& xv = mesh_->node_get_coordinates(v, &xv);
    xv -= xm;

    j = 2 * n;
    k = 2 * n + 1;

    a1(j) = 1.0;   a1(k) = 0.0;
    a2(j) = xv.y;  a2(k) =-xv.x;
    a3(j) = 0.0;   a3(k) = 1.0;

    a4(j) = xv.x;  a4(k) = 0.0;
    a5(j) = xv.y;  a5(k) = xv.x;
    a6(j) = 0.0;   a6(k) = xv.y;
  }

  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    xe = f->e_xm(i) - xm;
    l = 2 * nfaces + i;

    a1(l) = normal.x;
    a2(l) = p2D(xe.y, -xe.x) * normal;
    a3(l) = normal.y;

    a4(l) = xe.x * normal.x;
    a5(l) = p2D(xe.y, xe.x) * normal;
    a6(l) = xe.y * normal.y;
  }

  // calculate R K R^T / volume
  array2D<real> Mu(3,3);

  coefM(0, 0)               = (T0 % CT0) * area;
  coefM(0, 1) = coefM(1, 0) = (T0 % CT1) * area;
  coefM(0, 2) = coefM(2, 0) = (T0 % CT2) * area;

  coefM(1, 1)               = (T1 % CT1) * area; 
  coefM(1, 2) = coefM(2, 1) = (T1 % CT2) * area;
  coefM(2, 2)               = (T2 % CT2) * area; 

  coefM.inverse(3);

  for (i=0; i<nx; i++ ) { 
    v1(0) = b1(i);  
    v1(1) = b2(i);  
    v1(2) = b3(i); 

    for( j=i; j<nx; j++ ) {
      v2(0) = b1(j);  
      v2(1) = b2(j);  
      v2(2) = b3(j); 

      AmultB(coefM, v2, v3);
      localA(i, j) = v1 * v3;
    }
  }

  // calculate N
  b1 = 0;  
  b2 = 0;  
  b3 = 0;
  for (i=0; i<nv; i++ ) {
    j = (i + 1) % nv;

    lenght = f->e_lenght(i);
    normal = f->e_normal(i);
    tau = (f->xv(j) - f->xv(i)) / lenght;

    c = CT0 * normal;
    t = (tau    * c) * lenght / 2;
    p = (normal * c) * lenght;

    b1(2*i) += tau.x * t;  b1(2*i + 1) += tau.y * t;  b1(2*nv + i) += p;
    b1(2*j) += tau.x * t;  b1(2*j + 1) += tau.y * t;

    c = CT1 * normal;
    t = (tau    * c) * lenght / 2;
    p = (normal * c) * lenght;

    b2(2*i) += tau.x * t;  b2(2*i + 1) += tau.y * t;  b2(2*nv + i) += p;
    b2(2*j) += tau.x * t;  b2(2*j + 1) += tau.y * t;

    c = CT2 * normal;
    t = (tau    * c) * lenght / 2;
    p = (normal * c) * lenght;

    b3(2*i) += tau.x * t;  b3(2*i + 1) += tau.y * t;  b3(2*nv + i) += p;
    b3(2*j) += tau.x * t;  b3(2*j + 1) += tau.y * t;
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



