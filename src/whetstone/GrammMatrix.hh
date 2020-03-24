/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Calculation of various Gramm matrices. Implemented algorithms use
  integrator and optional database of monomial integrals.
*/

#ifndef AMANZI_WHETSTONE_GRAMM_MATRIX_HH_
#define AMANZI_WHETSTONE_GRAMM_MATRIX_HH_

#include <cmath>
#include <vector>

#include "Basis_Regularized.hh"
#include "Mesh.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"


namespace Amanzi {
namespace WhetStone {

class Polynomial;
class PolynomialOnMesh;

// Gramm matrix for polynomials
template<class MyMesh>
void GrammMatrix(
    NumericalIntegration<MyMesh>& numi,
    int order, PolynomialOnMesh& integrals,
    const Basis_Regularized<MyMesh>& basis, DenseMatrix& G);

// Gramm matrix for gradient of polynomials with tensorial weight
template<class MyMesh>
void GrammMatrixGradients(
    const Tensor& K, 
    NumericalIntegration<MyMesh>& numi,
    int order, PolynomialOnMesh& integrals,
    const Basis_Regularized<MyMesh>& basis, DenseMatrix& G);


/* ******************************************************************
* Gramm matrix G for polynomials.
****************************************************************** */
template<class MyMesh>
void GrammMatrix(
    NumericalIntegration<MyMesh>& numi,
    int order, PolynomialOnMesh& integrals,
    const Basis_Regularized<MyMesh>& basis, DenseMatrix& G)
{
  int d = numi.dimension();

  PolynomialIterator it0(d), it1(d);
  it0.begin(0);
  it1.begin(order + 1);

  int nd = it1.PolynomialPosition();
  G.Reshape(nd, nd);

  // extended database of integrals of monomials
  int c = integrals.id();
  numi.UpdateMonomialIntegralsCell(c, 2 * order, integrals);

  int multi_index[3];
  for (auto it = it0; it < it1; ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();
    double scalek = basis.monomial_scales()[it.MonomialSetOrder()];

    for (auto jt = it; jt < it1; ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      double scalel = basis.monomial_scales()[jt.MonomialSetOrder()];
      
      for (int i = 0; i < d; ++i) {
        multi_index[i] = index[i] + jndex[i];
      }

      int pos = PolynomialPosition(d, multi_index); 
      G(k, l) = G(l, k) = integrals.poly()(pos) * scalek * scalel; 
    }
  }
}


/* ******************************************************************
* Gramm matrix for gradient of polynomials with tensorial weight.
****************************************************************** */
template<class MyMesh>
void GrammMatrixGradients(
    const Tensor& K, 
    NumericalIntegration<MyMesh>& numi,
    int order, PolynomialOnMesh& integrals,
    const Basis_Regularized<MyMesh>& basis, DenseMatrix& G)
{
  int d = numi.dimension();

  PolynomialIterator it0(d), it1(d);
  it0.begin(0);
  it1.begin(order + 1);

  int nd = it1.PolynomialPosition();
  G.Reshape(nd, nd);

  // extended database of integrals of monomials
  int c = integrals.id();
  numi.UpdateMonomialIntegralsCell(c, std::max(0, 2 * (order - 1)), integrals);

  Tensor Ktmp(d, 2);
  if (K.rank() == 2)
    Ktmp = K;
  else 
    Ktmp.MakeDiagonal(K(0, 0));

  int multi_index[3];
  for (auto it = it0; it < it1; ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();
    double scalek = basis.monomial_scales()[it.MonomialSetOrder()];

    for (auto jt = it; jt < it1; ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      double scalel = basis.monomial_scales()[jt.MonomialSetOrder()];
      
      for (int i = 0; i < d; ++i) {
        multi_index[i] = index[i] + jndex[i];
      }

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

}  // namespace WhetStone
}  // namespace Amanzi

#endif

