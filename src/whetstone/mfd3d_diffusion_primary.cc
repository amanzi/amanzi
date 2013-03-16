/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: ara-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#include <cmath>
#include <vector>

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for stiffness matrix in heat conduction. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D::H1consistency(int cell, const Tensor& T,
                         Teuchos::SerialDenseMatrix<int, double>& N,
                         Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(cell, &nodes);
  int num_nodes = nodes.size();
  if (num_nodes != N.numRows()) return num_nodes;  // matrix was not reshaped

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  /* to calculate matrix R, we use temporary matrix N */
  N = 0;

  int num_faces = faces.size();
  for (int i = 0; i < num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    AmanziMesh::Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int num_face_nodes = face_nodes.size();

    for (int j = 0; j < num_face_nodes; j++) {
      int jnext = (j + 1) % num_face_nodes;
      int jprev = (j + num_face_nodes - 1) % num_face_nodes;

      int v = face_nodes[j];
      int vnext = face_nodes[jnext];
      int vprev = face_nodes[jprev];

      mesh_->node_get_coordinates(v, &p);
      mesh_->node_get_coordinates(vnext, &pnext);
      mesh_->node_get_coordinates(vprev, &pprev);

      v1 = pprev - pnext;
      v2 = p - fm;
      v3 = v1^v2;
      double u = dirs[i] * norm(v3) / (4 * area);

      int pos = FindPosition_(v, nodes);
      for (int k = 0; k < d; k++) N(pos, k) += normal[k] * u;
    }
  }

  for (int i = 0; i < num_nodes; i++) {  // calculate R T R^T / volume
    for (int k = 0; k < d; k++) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j = i; j < num_nodes; j++) {
      for (int k = 0; k < d; k++) v1[k] = N(j, k);
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);
  for (int i = 0; i < num_nodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    for (int k = 0; k < d; k++) N(i, k) = p[k] - cm[k];
    N(i, d) = 1;  // additional colum is added to the consistency condition
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



