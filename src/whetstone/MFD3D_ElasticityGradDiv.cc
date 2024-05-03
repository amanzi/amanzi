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
#include "MFD3D_ElasticityGradDiv.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Stress recostruction from nodal values
****************************************************************** */
int
MFD3D_ElasticityGradDiv::StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
{
  const auto& nodes = mesh_->getCellNodes(c);
  int nnodes = nodes.size();

  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  int nrows = d_ * nnodes;
  DenseVector div(nrows);
  div.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const auto& fnodes = mesh_->getFaceNodes(f);
    int mnodes = fnodes.size();

    const auto& normal = mesh_->getFaceNormal(f);
    double a = mesh_->getFaceArea(f) * dirs[n] / mnodes;

    for (int m = 0; m < mnodes; ++m) {
      int v = fnodes(m);
      int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));

      for (int k = 0; k < d_; ++k) {
        div(d_ * pos + k) += a * normal[k];
      }
    }
  }

  A.Reshape(nrows, nrows);

  double coef = T(0, 0);
  for (int m = 0; m < nrows; ++m) {
    for (int n = m; n < nrows; ++n) {
      A(n, m) = A(m, n) = coef * div(m) * div(n);
    }
  }

  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
