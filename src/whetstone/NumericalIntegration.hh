/*
  WhetStone, Version 2.2
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

#include <memory>

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "Monomial.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Quadrature1D.hh"
#include "Quadrature2D.hh"
#include "Quadrature3D.hh"
#include "SurfaceCoordinateSystem.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

// non-member functions
double
IntegrateFunctionsTriangle(const std::vector<AmanziGeometry::Point>& xy,
                           const std::vector<const WhetStoneFunction*>& funcs,
                           int order);

double
IntegrateFunctionsTetrahedron(const std::vector<AmanziGeometry::Point>& xy,
                              const std::vector<const WhetStoneFunction*>& funcs,
                              int order);

double
IntegrateFunctionsEdge(const AmanziGeometry::Point& x1,
                       const AmanziGeometry::Point& x2,
                       const std::vector<const WhetStoneFunction*>& funcs,
                       int order);

double
IntegratePolynomialsEdge(const AmanziGeometry::Point& x1,
                         const AmanziGeometry::Point& x2,
                         const std::vector<const PolynomialBase*>& polys);


class NumericalIntegration {
 public:
  NumericalIntegration(Teuchos::RCP<const AmanziMesh::MeshLight> mesh);
  ~NumericalIntegration(){};

  // main methods
  // integrate product of functions with quadrature order
  double IntegrateFunctionsTriangulatedCell(int c,
                                            const std::vector<const WhetStoneFunction*>& funcs,
                                            int order) const;

  double IntegrateFunctionsTriangulatedFace(int f,
                                            const std::vector<const WhetStoneFunction*>& funcs,
                                            int order) const;

  double IntegrateFunctionsEdge(const AmanziGeometry::Point& x1,
                                const AmanziGeometry::Point& x2,
                                const std::vector<const WhetStoneFunction*>& funcs,
                                int order) const;

  double
  IntegrateFunctionsEdge(int e, const std::vector<const WhetStoneFunction*>& funcs, int order) const
  {
    int v1, v2;
    AmanziGeometry::Point x1(d_), x2(d_);
    mesh_->edge_get_nodes(e, &v1, &v2);
    mesh_->node_get_coordinates(v1, &x1);
    mesh_->node_get_coordinates(v2, &x2);
    return WhetStone::IntegrateFunctionsEdge(x1, x2, funcs, order);
  }

  // integrate product of polynomials and monomials with different origins
  double IntegratePolynomialsCell(int c, const std::vector<const PolynomialBase*>& polys);

  double IntegratePolynomialsCell(int c,
                                  const std::vector<const PolynomialBase*>& polys,
                                  PolynomialOnMesh& integrals) const;

  double IntegratePolynomialsFace(int f, const std::vector<const PolynomialBase*>& polys) const;

  double IntegratePolynomialsEdge(int e, const std::vector<const PolynomialBase*>& polys) const
  {
    int v1, v2;
    AmanziGeometry::Point x1(d_), x2(d_);
    mesh_->edge_get_nodes(e, &v1, &v2);
    mesh_->node_get_coordinates(v1, &x1);
    mesh_->node_get_coordinates(v2, &x2);
    return WhetStone::IntegratePolynomialsEdge(x1, x2, polys);
  }

  // integral over a simplex
  double IntegrateFunctionsSimplex(const std::vector<AmanziGeometry::Point>& xy,
                                   const std::vector<const WhetStoneFunction*>& funcs,
                                   int order) const
  {
    if (xy.size() == 3)
      return IntegrateFunctionsTriangle(xy, funcs, order);
    else if (xy.size() == 4)
      return IntegrateFunctionsTetrahedron(xy, funcs, order);
    return 0.0;
  }

  // integrate group of monomials
  void IntegrateMonomialsCell(int c, int k, Polynomial& integrals) const;
  void UpdateMonomialIntegralsCell(int c, int order, Polynomial& integrals) const;
  void UpdateMonomialIntegralsCell(int c, int order, PolynomialOnMesh& integrals) const;

  // useful functions: integrate single polynomial
  double IntegratePolynomialCell(int c, const Polynomial& poly);

  double IntegratePolynomialFace(int f, const Polynomial& poly) const
  {
    const std::vector<const PolynomialBase*> polys(1, &poly);
    return IntegratePolynomialsFace(f, polys);
  }

  double IntegratePolynomialEdge(const AmanziGeometry::Point& x1,
                                 const AmanziGeometry::Point& x2,
                                 const Polynomial& poly) const
  {
    const std::vector<const PolynomialBase*> polys(1, &poly);
    return WhetStone::IntegratePolynomialsEdge(x1, x2, polys);
  }

  // useful functions: integrate single function
  double IntegrateFunctionEdge(const AmanziGeometry::Point& x1,
                               const AmanziGeometry::Point& x2,
                               const WhetStoneFunction& func,
                               int order) const
  {
    const std::vector<const WhetStoneFunction*> funcs(1, &func);
    return IntegrateFunctionsEdge(x1, x2, funcs, order);
  }

  // various bounds
  double PolynomialMaxValue(int f, const Polynomial& poly);

 private:
  void IntegrateMonomialsFace_(int c, int f, double factor, int k, Polynomial& integrals) const;

  void
  IntegrateMonomialsFaceReduction_(int c, int f, double factor, int k, Polynomial& integrals) const;

  void IntegrateMonomialsEdge_(const AmanziGeometry::Point& x1,
                               const AmanziGeometry::Point& x2,
                               double factor,
                               int k,
                               Polynomial& integrals) const;

 private:
  Teuchos::RCP<const AmanziMesh::MeshLight> mesh_;
  int d_;

  // often, we need to perform multiple integrations on the same object.
  // saving database of monomial integrals reduces computational cost,
  // especially in 3D
  PolynomialOnMesh integrals_;
};


//  -- trace of polynomials on manifold
Polynomial
ConvertPolynomialsToSurfacePolynomial(const AmanziGeometry::Point& xf,
                                      const std::shared_ptr<SurfaceCoordinateSystem>& coordsys,
                                      const std::vector<const PolynomialBase*>& polys);

} // namespace WhetStone
} // namespace Amanzi

#endif
