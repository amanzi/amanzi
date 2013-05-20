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

#include "mfd3d_elasticity.hh"
#include "tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for stifness matrix in geomechanics. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_Elasticity::H1consistency(int cell, const Tensor& T,
                                    Teuchos::SerialDenseMatrix<int, double>& N,
                                    Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  int nrows = N.numRows();

  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(cell, &nodes);
  int num_nodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int d = mesh_->space_dimension();
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  // convolution of tensors
  std::vector<Tensor> TE;

  for (int k = 0; k < d; k++) {
    Tensor E(d, 2);
    E(k, k) = 1.0;
    TE.push_back(T * E);
  }

  for (int k = 0; k < d; k++) {
    for (int l = k + 1; l < d; l++) {
      Tensor E(d, 2);
      E(k, l) = E(l, k) = 1.0;
      TE.push_back(T * E);
    }
  }
  int nd = TE.size();

  // to calculate matrix R, we use temporary matrix N
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
      int v = face_nodes[j];
      double u(0.5);

      if (d == 2) {
        u = 0.5 * dirs[i]; 
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes[jnext];
        int vprev = face_nodes[jprev];

        mesh_->node_get_coordinates(v, &p);
        mesh_->node_get_coordinates(vnext, &pnext);
        mesh_->node_get_coordinates(vprev, &pprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1^v2;
        u = dirs[i] * norm(v3) / (4 * area);
      }

      int pos = FindPosition_(v, nodes);
      for (int k = 0; k < nd; k++) {
        v1 = TE[k] * normal;
        for (int l = 0; l < d; l++) N(l * num_nodes + pos, k) += v1[l] * u;
      }
    }
  }

  // calculate R inv(T) R^T / volume
  Tensor Tinv(T);
  Tinv.inverse();

  double volume = mesh_->cell_volume(cell);
  Tinv *= 1.0 / volume;

  Teuchos::SerialDenseMatrix<int, double> RT(nrows, nd);  // R = N at this point
  if (Tinv.rank() == 1) {
    double* data_N = N.values();
    double* data_RT = RT.values();
    double s = Tinv(0, 0);
    for (int i = 0; i < nrows * nd; i++) data_RT[i] = data_N[i] * s;
  } else if (Tinv.rank() == 4) {
    Teuchos::SerialDenseMatrix<int, double> Ttmp(Teuchos::View, Tinv.data(), nd, nd, nd);
    MatrixMatrixProduct_(N, Ttmp, false, RT);
  }
  MatrixMatrixProduct_(RT, N, true, Ac); 

  // calculate matrix N
  N = 0.0;
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  for (int i = 0; i < num_nodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    v1 = p - cm;

    int md = 0;
    for (int k = 0; k < d; k++) {
      N(k * num_nodes + i, md) = v1[k];
      md++;
    }
    for (int k = 0; k < d; k++) {
      for (int l = k + 1; l < d; l++) {
        N(k * num_nodes + i, md) = v1[l];
        N(l * num_nodes + i, md) = v1[k];
        md++;
      }
    }
    for (int k = 0; k < d; k++) {  // additional columns correspod to kernel
      N(k * num_nodes + i, md) = 1.0;
      md++;
    }
    for (int k = 0; k < d; k++) {
      for (int l = k + 1; l < d; l++) {
        N(k * num_nodes + i, md) =  v1[l];
        N(l * num_nodes + i, md) = -v1[k];
        md++;
      }
    }
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrix(int cell, const Tensor& deformation,
                                      Teuchos::SerialDenseMatrix<int, double>& A)
{
  int d = mesh_->space_dimension();
  int nd = d * (d + 1);
  int nrows = A.numRows();

  Teuchos::SerialDenseMatrix<int, double> N(nrows, nd);
  Teuchos::SerialDenseMatrix<int, double> Ac(nrows, nrows);

  int ok = H1consistency(cell, deformation, N, Ac);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(cell, N, Ac, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Classical matrix-matrix product
****************************************************************** */
void MFD3D_Elasticity::MatrixMatrixProduct_(
    const Teuchos::SerialDenseMatrix<int, double>& A,
    const Teuchos::SerialDenseMatrix<int, double>& B, bool transposeB,
    Teuchos::SerialDenseMatrix<int, double>& AB)
{
  int nrows = A.numRows();
  int ncols = A.numCols();

  int mrows = B.numRows();
  int mcols = B.numCols();

  if (transposeB) {
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < mrows; j++) {
        double s = 0.0;
        for (int k = 0; k < ncols; k++) s += A(i, k) * B(j, k);
        AB(i, j) = s;
      }
    }
  } else {
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < mcols; j++) {
        double s = 0.0;
        for (int k = 0; k < ncols; k++) s += A(i, k) * B(k, j);
        AB(i, j) = s;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi



