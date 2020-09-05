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

#include <memory>

#include "Teuchos_RCP.hpp"

#include "NumericalIntegration.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Integrate a product of functions over a 2D or 3D triangle.
****************************************************************** */
double IntegrateFunctionsTriangle(
    const std::vector<AmanziGeometry::Point>& xy,
    const std::vector<const WhetStoneFunction*>& funcs, int order)
{
  // calculate minimal quadrature rule 
  int m(order);
  AMANZI_ASSERT(m < 14);

  int n1 = q2d_order[m][1];
  int n2 = n1 + q2d_order[m][0];

  AmanziGeometry::Point y1 = xy[1] - xy[0];
  AmanziGeometry::Point y2 = xy[2] - xy[0];

  double integral(0.0);
  for (int n = n1; n < n2; ++n) { 
    auto ym = xy[0] + y1 * q2d_points[n][1] + y2 * q2d_points[n][2];

    double a(q2d_weights[n]);
    for (int i = 0; i < funcs.size(); ++i) {
      a *= funcs[i]->Value(ym, 0.0);
    }
    integral += a;      
  }

  double area = norm(y1^y2) / 2;

  return integral * area;
}


/* ******************************************************************
* Integrate a product of functions over a tetrahedron
****************************************************************** */
double IntegrateFunctionsTetrahedron(
    const std::vector<AmanziGeometry::Point>& xy,
    const std::vector<const WhetStoneFunction*>& funcs, int order)
{
  int d(xy[0].dim());

  // calculate minimal quadrature rule 
  int m(order);
  AMANZI_ASSERT(m < 7);

  int n1 = q3d_order[m][1];
  int n2 = n1 + q3d_order[m][0];

  AmanziGeometry::Point ym(d);
  AmanziGeometry::Point y1 = xy[1] - xy[0];
  AmanziGeometry::Point y2 = xy[2] - xy[0];
  AmanziGeometry::Point y3 = xy[3] - xy[0];

  double integral(0.0);
  for (int n = n1; n < n2; ++n) { 
    ym = xy[0] + y1 * q3d_points[n][1] + y2 * q3d_points[n][2] + y3 * q3d_points[n][3];

    double a(q3d_weights[n]);
    for (int i = 0; i < funcs.size(); ++i) {
      a *= funcs[i]->Value(ym, 0.0);
    }
    integral += a;      
  }

  double volume = std::fabs(((y1^y2) * y3) / 6);

  return integral * volume;
}


/* ******************************************************************
* Integrate over face f a group of non-normalized monomials of
* the same order k centered at the centroid of cell c.
****************************************************************** */
template <>
void NumericalIntegration<AmanziMesh::Mesh>::IntegrateMonomialsFaceReduction_(
    int c, int f, double factor, int k, Polynomial& integrals) const
{
  int nk = PolynomialSpaceDimension(d_, k - 1);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);

  // create a surface mesh
  auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);
  Teuchos::RCP<const SurfaceMiniMesh> surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));
  NumericalIntegration<SurfaceMiniMesh> numi_f(surf_mesh);

  PolynomialIterator it(d_);
  for (it.begin(k); it.MonomialSetOrder() <= k; ++it) {
    int l = it.MonomialSetPosition();

    // using monomial centered at xc, create 2D polynomial centered at xf
    const int* idx = it.multi_index();
    Polynomial poly(d_, idx, 1.0);
    poly.set_origin(xc);
    poly.ChangeCoordinates(xf, *coordsys->tau());  

    integrals(nk + l) += factor * numi_f.IntegratePolynomialCell(f, poly);
  }
}


/* ******************************************************************
* Integrate over edge (x1,x2) a product of functions.
****************************************************************** */
double IntegrateFunctionsEdge(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<const WhetStoneFunction*>& funcs, int order)
{
  int m = order / 2;
  AMANZI_ASSERT(m < 8);

  int d = x1.dim();
  AmanziGeometry::Point xm(d);

  double integral(0.0);
  for (int n = 0; n <= m; ++n) { 
    double q1d(q1d_points[m][n]);
    for (int i = 0; i < d; ++i) {
      xm[i] = x1[i] * q1d + x2[i] * (1.0 - q1d);
    }

    double a(q1d_weights[m][n]);
    for (int i = 0; i < funcs.size(); ++i) {
      a *= funcs[i]->Value(xm, 0.0);
    }
    integral += a;      
  }

  return integral * norm(x2 - x1);
}


/* ******************************************************************
* Integrate over edge (x1,x2) a product of polynomials that may have
* different origins.
****************************************************************** */
double IntegratePolynomialsEdge(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<const PolynomialBase*>& polys)
{
  // minimal quadrature rule
  int k(0);
  for (int i = 0; i < polys.size(); ++i) {
    k += polys[i]->order();
  }
  int m = k / 2;
  AMANZI_ASSERT(m < 8);

  AmanziGeometry::Point xm(x1.dim());

  double integral(0.0);
  for (int n = 0; n <= m; ++n) { 
    xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);

    double a(q1d_weights[m][n]);
    for (int i = 0; i < polys.size(); ++i) {
      a *= polys[i]->Value(xm);
    }
    integral += a;      
  }

  return integral * norm(x2 - x1);
}


/* ******************************************************************
* TBW
****************************************************************** */
template <>
void NumericalIntegration<SurfaceMiniMesh>::IntegrateMonomialsFaceReduction_(
    int c, int f, double factor, int k, Polynomial& integrals) const
{};


/* ******************************************************************
* Create surface integration object
****************************************************************** */
Polynomial ConvertPolynomialsToSurfacePolynomial(
    const AmanziGeometry::Point& xf, 
    const std::shared_ptr<SurfaceCoordinateSystem>& coordsys,
    const std::vector<const PolynomialBase*>& polys)
{
  int d = xf.dim();
  Polynomial product(d - 1, 0);
  product(0) = 1.0;

  for (int i = 0; i < polys.size(); ++ i) {
    Polynomial tmp(d, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->get_origin());
    tmp.ChangeCoordinates(xf, *coordsys->tau());  
    product *= tmp;
  }

  return product;
}

}  // namespace WhetStone
}  // namespace Amanzi

