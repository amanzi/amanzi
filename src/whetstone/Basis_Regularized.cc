/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The regularized basis for dG methods: x^k y^l / h^(k+l), where
  h is a measure of cell size.
*/

#include "Basis_Regularized.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Prepare scaling data for the regularized basis.
****************************************************************** */
void Basis_Regularized::Init(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int c, int order)
{
  int k0 = monomial_scales_.size();

  if (k0 < order + 1) {
    int d = mesh->space_dimension();
    double volume = mesh->cell_volume(c);

    monomial_scales_.resize(order + 1);
    for (int k = k0; k < order + 1; ++k) {
      monomial_scales_[k] = std::pow(volume, -(double)k / d);
    }
  }
}


/* ******************************************************************
* Recover polynomial from data coeffieints. 
****************************************************************** */
Polynomial Basis_Regularized::CalculatePolynomial(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    int c, int order, DenseVector& coefs) const
{
  int d = mesh->space_dimension();
  Polynomial poly(d, order);

  poly.SetPolynomialCoefficients(coefs);
  poly.set_origin(mesh->cell_centroid(c));

  for (int k = 0; k < order + 1; ++k) {
    auto& mono = poly.MonomialSet(k);
    mono *= monomial_scales_[k];
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

