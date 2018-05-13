/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Numerical and exact integration over polytopal cells. An object of
  this class can be used both for the whole mesh and single cell. 
*/

#ifndef AMANZI_WHETSTONE_NUMERICAL_INTEGRATION_HH_
#define AMANZI_WHETSTONE_NUMERICAL_INTEGRATION_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Natural.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Quadrature1D.hh"
#include "Quadrature2D.hh"

namespace Amanzi {
namespace WhetStone {

class NumericalIntegration { 
 public:
  NumericalIntegration(Teuchos::RCP<const AmanziMesh::Mesh> mesh, bool single_cell);
  ~NumericalIntegration() {};

  // main methods
  // integrate product of polynomials with potentialy different origins
  // -- automatically calculate quadrature order if order < 0 
  double IntegratePolynomialsTrianglatedCell(
      int c, const std::vector<const Polynomial*>& polys, int order = -1) const;

  double IntegratePolynomialsFace(
      int f, const std::vector<const Polynomial*>& polys) const;

  double IntegratePolynomialsEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const std::vector<const Polynomial*>& polys) const;

  // -- automatically calculate quadrature order if order < 0
  double IntegratePolynomialsTriangle(
      const std::vector<AmanziGeometry::Point>& xy,
      const std::vector<const Polynomial*>& polys, int order = -1) const;

  double IntegratePolynomialsTriangle(
      const std::vector<AmanziGeometry::Point>& xy, const Polynomial& poly) const {
    const std::vector<const Polynomial*> polys(1, &poly);
    return IntegratePolynomialsTriangle(xy, polys);
  }

  // integrate group of monomials 
  void IntegrateMonomialsCell(int c, MonomialSet& monomials);
  void UpdateMonomialIntegralsCell(int c, int order, PolynomialOnMesh& integrals);

  // useful functions: integrate single polynomial
  double IntegratePolynomialCell(int c, const Polynomial& poly);

  double IntegratePolynomialFace(int f, const Polynomial& poly) const {
    const std::vector<const Polynomial*> polys(1, &poly);
    return IntegratePolynomialsFace(f, polys);
  }

  double IntegratePolynomialEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const Polynomial& poly) const {
    const std::vector<const Polynomial*> polys(1, &poly);
    return IntegratePolynomialsEdge(x1, x2, polys);
  }

  // various bounds
  double PolynomialMaxValue(int f, const Polynomial& poly);

  // natural scaling of monomials (e.g. x^k / h^k)
  // -- scaling factor is constant for monomials of the same order
  double MonomialNaturalScales(int c, int k);
  // -- polynomial is converted from regular to natural basis
  void ChangeBasisRegularToNatural(int c, Polynomial& p);
  // -- polynomial is converted from natural to regular basis
  void ChangeBasisNaturalToRegular(int c, Polynomial& p);

 private:
  void IntegrateMonomialsFace_(int c, int f, double factor, MonomialSet& monomials);
  void IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      double factor, MonomialSet& monomials);

  double MonomialNaturalSingleScale_(int k, double volume) const;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_;

  // cached variables
  bool single_cell_;
  Basis_Natural single_cell_basis_;
  std::vector<Basis_Natural> basis_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

