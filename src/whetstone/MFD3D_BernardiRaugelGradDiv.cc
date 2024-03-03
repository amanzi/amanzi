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

  The mimetic finite difference method for elasticity.
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_BernardiRaugelGradDiv.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem>
MFD3D_BernardiRaugelGradDiv::schema() const
{
  std::vector<SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::NODE, DOF_Type::POINT, d_));
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::FACE, DOF_Type::NORMAL_COMPONENT, 1));
  return items;
}


/* ******************************************************************
* Stress recostruction from nodal values
****************************************************************** */
int
MFD3D_BernardiRaugelGradDiv::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  const auto& nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  DenseVector div(nfaces);
  div.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    div(n) = mesh_->getFaceArea(f) * dirs[n];
  }

  int nrows0 = d_ * nnodes;
  int nrows = nrows0 + nfaces;
  A.Reshape(nrows, nrows);

  double coef = T(0, 0);
  for (int m = 0; m < nfaces; ++m) {
    int m0 = nrows0 + m;
    for (int n = m; n < nfaces; ++n) {
      int n0 = nrows0 + n;
      A(n0, m0) = A(m0, n0) = coef * div(m) * div(n);
    }
  }

  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
