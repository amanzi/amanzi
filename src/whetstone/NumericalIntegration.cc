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

#include <cmath>

#include "Monomial.hh"
#include "NumericalIntegration.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Constructor.
 ****************************************************************** */
NumericalIntegration::NumericalIntegration(
  Teuchos::RCP<const AmanziMesh::Mesh> mesh)
  : mesh_(mesh), d_(mesh->space_dimension()){};


/* ******************************************************************
 * Integrate over triangulated cell a product of functions.
 ****************************************************************** */
double
NumericalIntegration::IntegrateFunctionsTriangulatedCell(
  int c, const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  double integral(0.0);

  Kokkos::View<AmanziMesh::Entity_ID*> faces, nodes;
  std::vector<AmanziGeometry::Point> xy(d_ + 1);

  mesh_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);

  xy[0] = mesh_->cell_centroid(c);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces(n);
    mesh_->face_get_nodes(f, nodes);
    int nnodes = nodes.extent(0);

    if (d_ == 3) {
      xy[1] = mesh_->face_centroid(f);

      for (int k = 0; k < nnodes; ++k) {
        int l = (k + 1) % nnodes;
        mesh_->node_get_coordinates(nodes(k), &(xy[2]));
        mesh_->node_get_coordinates(nodes(l), &(xy[3]));

        integral += IntegrateFunctionsTetrahedron_(xy, funcs, order);
      }
    } else if (d_ == 2) {
      mesh_->node_get_coordinates(nodes(0), &(xy[1]));
      mesh_->node_get_coordinates(nodes(1), &(xy[2]));

      integral += IntegrateFunctionsTriangle_(xy, funcs, order);
    }
  }

  return integral;
}


/* ******************************************************************
 * Integrate over triangulated face a product of functions.
 ****************************************************************** */
double
NumericalIntegration::IntegrateFunctionsTriangulatedFace(
  int f, const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  double integral(0.0);
  AmanziGeometry::Point x1(d_), x2(d_);
  Kokkos::View<AmanziMesh::Entity_ID*> faces, nodes;

  if (d_ == 3) {
    std::vector<AmanziGeometry::Point> xy(d_ + 1);

    mesh_->face_get_nodes(f, nodes);
    int nnodes = nodes.extent(0);

    xy[0] = mesh_->face_centroid(f);

    for (int k = 0; k < nnodes; ++k) {
      int l = (k + 1) % nnodes;
      mesh_->node_get_coordinates(nodes(k), &(xy[1]));
      mesh_->node_get_coordinates(nodes(l), &(xy[2]));

      integral += IntegrateFunctionsTriangle_(xy, funcs, order);
    }
  } else if (d_ == 2) {
    mesh_->face_get_nodes(f, nodes);

    mesh_->node_get_coordinates(nodes(0), &x1);
    mesh_->node_get_coordinates(nodes(1), &x2);

    integral += IntegrateFunctionsEdge(x1, x2, funcs, order);
  }

  return integral;
}


/* ******************************************************************
 * Integrate over edge (x1,x2) a product of functions.
 ****************************************************************** */
double
NumericalIntegration::IntegrateFunctionsEdge(
  const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
  const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  int m = order / 2;
  AMANZI_ASSERT(m < 8);

  AmanziGeometry::Point xm(d_);

  double integral(0.0);
  for (int n = 0; n <= m; ++n) {
    double q1d(q1d_points[m][n]);
    for (int i = 0; i < d_; ++i) { xm[i] = x1[i] * q1d + x2[i] * (1.0 - q1d); }

    double a(q1d_weights[m][n]);
    for (int i = 0; i < funcs.size(); ++i) { a *= funcs[i]->Value(xm); }
    integral += a;
  }

  return integral * norm(x2 - x1);
}


/* ******************************************************************
 * Integrate polynomial over cell c.
 ****************************************************************** */
double
NumericalIntegration::IntegratePolynomialCell(int c, const Polynomial& poly)
{
  // calculate integrals of monomials centered at cell centroid
  int order = poly.order();
  Polynomial integrals(d_, order);

  for (int k = 0; k <= order; ++k) { IntegrateMonomialsCell(c, k, integrals); }

  // dot product of coefficients of two polynomials.
  Polynomial tmp = poly;
  tmp.ChangeOrigin(mesh_->cell_centroid(c));

  double value(0.0);
  for (int n = 0; n < tmp.size(); ++n) { value += integrals(n) * tmp(n); }

  return value;
}


/* ******************************************************************
 * Integrate product of polynomials and monomials over cells c. They
 * may have different origins.
 ****************************************************************** */
double
NumericalIntegration::IntegratePolynomialsCell(
  int c, const std::vector<const PolynomialBase*>& polys) const
{
  // create a single polynomial centered at cell centroid
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  Polynomial product(d_, 0);
  product(0) = 1.0;
  product.set_origin(xc);

  for (int i = 0; i < polys.size(); ++i) {
    Polynomial tmp(d_, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->origin());
    tmp.ChangeOrigin(xc);
    product *= tmp;
  }

  // calculate integrals of monomials centered at cell centroid
  int order = product.order();
  Polynomial integrals(d_, order);

  for (int k = 0; k <= order; ++k) { IntegrateMonomialsCell(c, k, integrals); }

  // dot product of coefficients of two polynomials.
  double value(0.0);
  for (int n = 0; n < product.size(); ++n) {
    value += integrals(n) * product(n);
  }

  return value;
}


/* ******************************************************************
 * Integrate over face f the product of polynomials and monomials that
 * may have different origins.
 ****************************************************************** */
double
NumericalIntegration::IntegratePolynomialsFace(
  int f, const std::vector<const PolynomialBase*>& polys) const
{
  AmanziGeometry::Point enormal(d_), x1(d_), x2(d_);

  if (d_ == 2) {
    Kokkos::View<Entity_ID*> nodes;
    mesh_->face_get_nodes(f, nodes);

    mesh_->node_get_coordinates(nodes(0), &x1);
    mesh_->node_get_coordinates(nodes(1), &x2);
    return IntegratePolynomialsEdge(x1, x2, polys);
  }

  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point fnormal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  fnormal /= area;

  // create a single polynomial centered at face centroid
  Polynomial product(d_, 0);
  product(0) = 1.0;
  product.set_origin(xf);

  for (int i = 0; i < polys.size(); ++i) {
    Polynomial tmp(d_, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->origin());
    tmp.ChangeOrigin(xf);
    product *= tmp;
  }

  // Apply Euler theorem to each monomial
  double sum(0.0);
  Kokkos::View<Entity_ID*> edges;
  Kokkos::View<int*> dirs;

  mesh_->face_get_edges_and_dirs(f, edges, &dirs);
  int nedges = edges.extent(0);

  for (int n = 0; n < nedges; ++n) {
    int e = edges(n);
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double length = mesh_->edge_length(e);

    enormal = tau ^ fnormal;

    // rescale polynomial coefficients
    double tmp = dirs(n) * ((xe - xf) * enormal) / length;

    Polynomial q(product);
    for (auto it = q.begin(); it < q.end(); ++it) {
      int m = it.MonomialSetOrder();
      int k = it.MonomialSetPosition();
      q(m, k) *= tmp / (m + 2);
    }

    // integrate along edge
    int n0, n1;
    mesh_->edge_get_nodes(e, &n0, &n1);
    mesh_->node_get_coordinates(n0, &x1);
    mesh_->node_get_coordinates(n1, &x2);

    std::vector<const PolynomialBase*> q_ptr(1, &q);
    sum += IntegratePolynomialsEdge(x1, x2, q_ptr);
  }

  return sum;
}


/* ******************************************************************
 * Integrate over edge (x1,x2) a product of polynomials that may have
 * different origins.
 ****************************************************************** */
double
NumericalIntegration::IntegratePolynomialsEdge(
  const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
  const std::vector<const PolynomialBase*>& polys) const
{
  // minimal quadrature rule
  int k(0);
  for (int i = 0; i < polys.size(); ++i) { k += polys[i]->order(); }
  int m = k / 2;
  AMANZI_ASSERT(m < 8);

  AmanziGeometry::Point xm(d_);

  double integral(0.0);
  for (int n = 0; n <= m; ++n) {
    xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);

    double a(q1d_weights[m][n]);
    for (int i = 0; i < polys.size(); ++i) { a *= polys[i]->Value(xm); }
    integral += a;
  }

  return integral * norm(x2 - x1);
}


/* ******************************************************************
 * Integrate a product of functions over a 2D or 3D triangle.
 ****************************************************************** */
double
NumericalIntegration::IntegrateFunctionsTriangle_(
  const std::vector<AmanziGeometry::Point>& xy,
  const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  // calculate minimal quadrature rule
  int m(order);
  AMANZI_ASSERT(m < 10);

  int n1 = q2d_order[m][1];
  int n2 = n1 + q2d_order[m][0];

  AmanziGeometry::Point y1 = xy[1] - xy[0];
  AmanziGeometry::Point y2 = xy[2] - xy[0];

  double integral(0.0);
  for (int n = n1; n < n2; ++n) {
    auto ym = xy[0] + y1 * q2d_points[n][1] + y2 * q2d_points[n][2];

    double a(q2d_weights[n]);
    for (int i = 0; i < funcs.size(); ++i) { a *= funcs[i]->Value(ym); }
    integral += a;
  }

  double area = norm(y1 ^ y2) / 2;

  return integral * area;
}


/* ******************************************************************
 * Integrate a product of functions over a tetrahedron
 ****************************************************************** */
double
NumericalIntegration::IntegrateFunctionsTetrahedron_(
  const std::vector<AmanziGeometry::Point>& xy,
  const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  // calculate minimal quadrature rule
  int m(order);
  AMANZI_ASSERT(m < 7);

  int n1 = q3d_order[m][1];
  int n2 = n1 + q3d_order[m][0];

  AmanziGeometry::Point ym(d_);
  AmanziGeometry::Point y1 = xy[1] - xy[0];
  AmanziGeometry::Point y2 = xy[2] - xy[0];
  AmanziGeometry::Point y3 = xy[3] - xy[0];

  double integral(0.0);
  for (int n = n1; n < n2; ++n) {
    ym = xy[0] + y1 * q3d_points[n][1] + y2 * q3d_points[n][2] +
         y3 * q3d_points[n][3];

    double a(q3d_weights[n]);
    for (int i = 0; i < funcs.size(); ++i) { a *= funcs[i]->Value(ym); }
    integral += a;
  }

  double volume = std::fabs(((y1 ^ y2) * y3) / 6);

  return integral * volume;
}


/* ******************************************************************
 * Integrate over cell c a group of non-normalized monomials of
 * the same order centered at the centroid of c.
 ****************************************************************** */
void
NumericalIntegration::UpdateMonomialIntegralsCell(int c, int order,
                                                  Polynomial& integrals)
{
  int k0 = integrals.order();

  if (k0 < order) {
    integrals.Reshape(d_, order);

    for (int k = k0 + 1; k <= order; ++k)
      IntegrateMonomialsCell(c, k, integrals);
  }
}


void
NumericalIntegration::UpdateMonomialIntegralsCell(int c, int order,
                                                  PolynomialOnMesh& integrals)
{
  Polynomial& poly = integrals.poly();
  int k0 = poly.order();

  if (integrals.kind() != (Entity_kind)WhetStone::CELL || integrals.id() != c) {
    integrals.set_kind((Entity_kind)WhetStone::CELL);
    integrals.set_id(c);
    k0 = -1;
  }

  // add additional integrals of monomials
  if (k0 < order) {
    poly.Reshape(d_, order);

    for (int k = k0 + 1; k <= order; ++k) IntegrateMonomialsCell(c, k, poly);
  }
}


/* ******************************************************************
 * Integrate over cell c a group of non-normalized monomials of
 * the same order centered at the centroid of c.
 ****************************************************************** */
void
NumericalIntegration::IntegrateMonomialsCell(int c, int k,
                                             Polynomial& integrals) const
{
  int nk = PolynomialSpaceDimension(d_, k - 1);
  int mk = MonomialSpaceDimension(d_, k);
  for (int i = 0; i < mk; ++i) { integrals(nk + i) = 0.0; }

  Kokkos::View<Entity_ID*> faces, nodes;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int nfaces = faces.extent(0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces(n);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double tmp = dirs(n) * ((xf - xc) * normal) / (k + d_);

    if (d_ == 3) {
      tmp /= mesh_->face_area(f);
      IntegrateMonomialsFace_(c, f, tmp, k, integrals);
    } else if (d_ == 2) {
      mesh_->face_get_nodes(f, nodes);

      AmanziGeometry::Point x1(d_), x2(d_);
      mesh_->node_get_coordinates(nodes(0), &x1);
      mesh_->node_get_coordinates(nodes(1), &x2);

      x1 -= xc; // simple change of origin
      x2 -= xc;
      IntegrateMonomialsEdge_(x1, x2, tmp, k, integrals);
    }
  }
}


/* ******************************************************************
 * Integrate over face f a group of non-normalized monomials of
 * the same order k centered at the centroid of cell c.
 ****************************************************************** */
void
NumericalIntegration::IntegrateMonomialsFace_(int c, int f, double factor,
                                              int k,
                                              Polynomial& integrals) const
{
  int nk = PolynomialSpaceDimension(d_, k - 1);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point normal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  normal /= area;

  AmanziGeometry::Point fnormal(d_), x1(d_), x2(d_);
  std::vector<const PolynomialBase*> polys(1);

  PolynomialIterator it(d_);
  for (it.begin(k); it.MonomialSetOrder() <= k; ++it) {
    int l = it.MonomialSetPosition();

    // using monomial centered at xc, create polynomial centred at xf
    const int* idx = it.multi_index();

    Monomial mono(d_, idx, 1.0);
    mono.set_origin(xc);
    Polynomial poly = integrals.ChangeOrigin(mono, xf);

    Kokkos::View<Entity_ID*> edges;
    Kokkos::View<int*> dirs;

    mesh_->face_get_edges_and_dirs(f, edges, &dirs);
    int nedges = edges.extent(0);

    for (int n = 0; n < nedges; ++n) {
      int e = edges(n);
      const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
      const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
      double length = mesh_->edge_length(e);

      fnormal = tau ^ normal;

      // rescale polynomial coefficients
      double tmp = (factor * dirs(n)) * ((xe - xf) * fnormal) / length;

      Polynomial q(poly);
      for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int s = jt.PolynomialPosition();
        q(s) *= tmp / (m + d_ - 1);
      }

      // integrate along edge
      int n0, n1;
      mesh_->edge_get_nodes(e, &n0, &n1);
      mesh_->node_get_coordinates(n0, &x1);
      mesh_->node_get_coordinates(n1, &x2);

      polys[0] = &q;
      integrals(nk + l) += IntegratePolynomialsEdge(x1, x2, polys);
    }
  }
}


/* ******************************************************************
 * Integrate over edge (x1,x2) a group of non-normalized monomials of
 * the same order k centered at zero.
 ****************************************************************** */
void
NumericalIntegration::IntegrateMonomialsEdge_(const AmanziGeometry::Point& x1,
                                              const AmanziGeometry::Point& x2,
                                              double factor, int k,
                                              Polynomial& integrals) const
{
  int nk = PolynomialSpaceDimension(d_, k - 1);

  // minimal quadrature rule
  int m = k / 2;
  AMANZI_ASSERT(m < 8);

  PolynomialIterator it(d_);
  for (it.begin(k); it.MonomialSetOrder() <= k; ++it) {
    const int* idx = it.multi_index();
    int l = it.MonomialSetPosition();

    for (int n = 0; n <= m; ++n) {
      double a1(factor), q1d(q1d_points[m][n]);
      for (int i = 0; i < d_; ++i) {
        double tmp = x1[i] * q1d + x2[i] * (1.0 - q1d);
        a1 *= std::pow(tmp, idx[i]);
      }

      integrals(nk + l) += a1 * q1d_weights[m][n];
    }
  }
}


/* ******************************************************************
 * Approximate maximum value of a polynomial.
 ****************************************************************** */
double
NumericalIntegration::PolynomialMaxValue(int f, const Polynomial& poly)
{
  int k = poly.order();
  int m = k / 2;
  double pmax;
  AmanziGeometry::Point x1(d_), x2(d_), xm(d_);

  Kokkos::View<Entity_ID*> nodes;
  mesh_->face_get_nodes(f, nodes);

  if (d_ == 2) {
    mesh_->node_get_coordinates(nodes(0), &x1);
    mesh_->node_get_coordinates(nodes(1), &x2);

    pmax = std::max(fabs(poly.Value(x1)), fabs(poly.Value(x2)));
    for (int n = 0; n <= m; ++n) {
      xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);
      pmax = std::max(pmax, fabs(poly.Value(xm)));
    }
  } else {
    pmax = -1.0e+99;
    for (int i = 0; i < nodes.extent(0); ++i) {
      mesh_->node_get_coordinates(nodes(i), &xm);
      pmax = std::max(pmax, fabs(poly.Value(xm)));
    }
  }

  return pmax;
}

} // namespace WhetStone
} // namespace Amanzi
