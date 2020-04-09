/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Supporting routines: Gramm matrices.
*/

#ifndef AMANZI_WHETSTONE_GRAMM_MATRIX_HH_
#define AMANZI_WHETSTONE_GRAMM_MATRIX_HH_

#include <cmath>
#include <vector>

#include "Basis_Regularized.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;
class PolynomialOnMesh;

// Gramm matrix for polynomials
void
GrammMatrix(const Polynomial& poly, const PolynomialOnMesh& integrals,
            const Basis_Regularized& basis, DenseMatrix& G);

// Gramm matrix for gradient of polynomials with tensorial weight
void
GrammMatrixGradients(const Tensor& K, const Polynomial& poly,
                     const PolynomialOnMesh& integrals,
                     const Basis_Regularized& basis, DenseMatrix& G);

} // namespace WhetStone
} // namespace Amanzi

#endif
