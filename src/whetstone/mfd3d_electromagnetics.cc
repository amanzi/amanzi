/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Release name: ara-to.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "tensor.hh"
#include "mfd3d_electromagnetics.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for the mass matrix in electromagnetics.
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D_Electromagnetics::L2consistency(int c, const Tensor& T,
                                          DenseMatrix& N, DenseMatrix& Mc)
{
  int nrows = N.NumRows();

  Entity_ID_List edges, nodes;

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  int d = mesh_->space_dimension();
  AmanziGeometry::Point p(d);

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);
  for (int i = 0; i < nedges; i++) Mc(i, i) = 1.0;

  // Matrix N are simply transformed tangents
  AmanziGeometry::Point p0(d), p1(d), tau(d), v2(d);

  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);

    v2 = T * tau;
    for (int k = 0; k < d; ++k) N(i, k) = v2[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Matrix matrix for edge-based discretization.
****************************************************************** */
int MFD3D_Electromagnetics::MassMatrix(int c, const Tensor& T, DenseMatrix& M)
{
  int d = mesh_->space_dimension();
  int nrows = M.NumRows();

  DenseMatrix N(nrows, d);
  DenseMatrix Mc(nrows, nrows);

  int ok = L2consistency(c, T, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  StabilityScalar(c, N, Mc, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi



