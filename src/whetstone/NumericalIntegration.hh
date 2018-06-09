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

#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Quadrature1D.hh"
#include "Quadrature2D.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class NumericalIntegration { 
 public:
  NumericalIntegration(Teuchos::RCP<const AmanziMesh::Mesh> mesh);
  ~NumericalIntegration() {};

  // main methods
  // integrate product of functions with quadrature order 
  double IntegrateFunctionsTrianglatedCell(
      int c, const std::vector<const WhetStoneFunction*>& funcs, int order) const;

  double IntegratePolynomialsFace(
      int f, const std::vector<const Polynomial*>& polys) const;

  double IntegratePolynomialsEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const std::vector<const Polynomial*>& polys) const;

  // -- automatically calculate quadrature order if order < 0
  double IntegrateFunctionsTriangle(
      const std::vector<AmanziGeometry::Point>& xy,
      const std::vector<const WhetStoneFunction*>& funcs, int order) const;

  double IntegratePolynomialTriangle(
      const std::vector<AmanziGeometry::Point>& xy, const Polynomial& poly) const {
    const std::vector<const WhetStoneFunction*> funcs(1, &poly);
    return IntegrateFunctionsTriangle(xy, funcs, poly.order());
  }

  // integrate group of monomials 
  void IntegrateMonomialsCell(int c, int k, Polynomial& integrals);
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

 private:
  void IntegrateMonomialsFace_(int c, int f, double factor, int k, Polynomial& integrals);
  void IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      double factor, int k, Polynomial& integrals);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

