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

#include "GrammMatrix.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Gramm matrix G for polynomials.
 ****************************************************************** */
void
GrammMatrix(const Polynomial& poly, const PolynomialOnMesh& integrals,
            const Basis_Regularized& basis, DenseMatrix& G)
{
  int nd = poly.size();
  int d = poly.dimension();
  G.Reshape(nd, nd);

  int multi_index[3];
  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();
    double scalek = basis.monomial_scales()[it.MonomialSetOrder()];

    for (auto jt = it; jt < poly.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      double scalel = basis.monomial_scales()[jt.MonomialSetOrder()];

      for (int i = 0; i < d; ++i) { multi_index[i] = index[i] + jndex[i]; }

      int pos = PolynomialPosition(d, multi_index);
      G(k, l) = G(l, k) = integrals.poly()(pos) * scalek * scalel;
    }
  }
}

/* ******************************************************************
 * Gramm matrix for gradient of polynomials with tensorial weight.
 ****************************************************************** */
void
GrammMatrixGradients(const Tensor& K, const Polynomial& poly,
                     const PolynomialOnMesh& integrals,
                     const Basis_Regularized& basis, DenseMatrix& G)
{
  int nd = poly.size();
  int d = poly.dimension();
  G.Reshape(nd, nd);

  Tensor Ktmp(d, 2);
  if (K.rank() == 2)
    Ktmp = K;
  else
    Ktmp.MakeDiagonal(K(0, 0));

  int multi_index[3];
  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();
    double scalek = basis.monomial_scales()[it.MonomialSetOrder()];

    for (auto jt = it; jt < poly.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      double scalel = basis.monomial_scales()[jt.MonomialSetOrder()];

      for (int i = 0; i < d; ++i) { multi_index[i] = index[i] + jndex[i]; }

      double sum(0.0), tmp;
      for (int i = 0; i < d; ++i) {
        if (index[i] > 0) {
          multi_index[i] -= 1;

          for (int j = 0; j < d; ++j) {
            if (jndex[j] > 0) {
              multi_index[j] -= 1;

              int n = PolynomialPosition(d, multi_index);
              tmp = integrals.poly()(n);
              sum += Ktmp(i, j) * tmp * index[i] * jndex[j];

              multi_index[j] += 1;
            }
          }
          multi_index[i] += 1;
        }
      }

      G(l, k) = G(k, l) = sum * scalek * scalel;
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
