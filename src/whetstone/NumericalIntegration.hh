/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Numerical and exact integration over polytopal cells.
*/

#ifndef AMANZI_WHETSTONE_NUMERICAL_INTEGRATION_HH_
#define AMANZI_WHETSTONE_NUMERICAL_INTEGRATION_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

// Gauss quadrature on interval (0,1)
const double q1d_weights[5][5] = {
    1.0, 0.0, 0.0, 0.0, 0.0,
    0.5, 0.5, 0.0, 0.0, 0.0,
    0.277777777777778, 0.444444444444444, 0.277777777777778, 0.0, 0.0,
    0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727, 0.0,
    0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095
};
const double q1d_points[5][5] = {
    0.5, 0.0, 0.0, 0.0, 0.0,
    0.211324865405187, 0.788675134594813, 0.0, 0.0, 0.0,
    0.112701665379258, 0.5, 0.887298334620742, 0.0, 0.0,
    0.0694318442029737, 0.330009478207572, 0.669990521792428, 0.930568155797026, 0.0,
    0.0469100770306680, 0.230765344947158, 0.0, 0.769234655052841, 0.953089922969332
};

class NumericalIntegration { 
 public:
  NumericalIntegration(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : mesh_(mesh),
      d_(mesh_->space_dimension()) {};

  ~NumericalIntegration() {};

  // main methods
  // integrate product of polynomials with potentialy different origins
  double IntegratePolynomialsFace(
      int f, const std::vector<const Polynomial*>& polys) const;

  double IntegratePolynomialsEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const std::vector<const Polynomial*>& polys) const;

  // integrate group of monomials 
  void IntegrateMonomialsCell(int c, Monomial& monomials);

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

  // miscalleneous: order is not used yet
  void set_order(int order) { order_ = order; }

 private:
  void IntegrateMonomialsFace_(int c, int f, double factor, Monomial& monomials);
  void IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      double factor, Monomial& monomials);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int order_, d_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

