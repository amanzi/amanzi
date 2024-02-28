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
  AMANZI_ASSERT(d_ == 2);

  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  int nrows = d_ * nfaces;
  DenseVector div(nrows);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const auto& normal = mesh_->getFaceNormal(f);
    double a = mesh_->getFaceArea(f) * dirs[n] / 2;

    int n1 = n;
    int n2 = (n + 1) % nfaces;

    for (int k = 0; k < d_; ++k) {
      div(d_ * n1 + k) += a * normal[k];
      div(d_ * n2 + k) += a * normal[k];
    }
  }

  A.Reshape(nrows, nrows);

  double coef = T(0, 0);
  for (int m = 0; m < nrows; ++m) {
    for (int n = m; n < nrows; ++n) { A(n, m) = A(m, n) = coef * div(m) * div(n); }
  }

  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
