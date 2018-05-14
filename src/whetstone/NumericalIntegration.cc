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

#include "NumericalIntegration.hh"
#include "Monomial.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor.
****************************************************************** */
NumericalIntegration::NumericalIntegration(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, bool single_cell)
  : mesh_(mesh),
    d_(mesh_->space_dimension()),
    single_cell_(single_cell)
{
  if (! single_cell_) {
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    basis_.resize(ncells_wghost);
  }
}


/* ******************************************************************
* Integrate over triangulated face a product of polynomials.
****************************************************************** */
double NumericalIntegration::IntegratePolynomialsTrianglatedCell(
    int c, const std::vector<const Polynomial*>& polys, int order) const
{
  double integral(0.0);

  if (d_ == 2) {
    AmanziMesh::Entity_ID_List faces, nodes;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    std::vector<AmanziGeometry::Point> xy(3); 
    xy[0] = mesh_->cell_centroid(c);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      mesh_->face_get_nodes(f, &nodes);
      mesh_->node_get_coordinates(nodes[0], &(xy[1]));
      mesh_->node_get_coordinates(nodes[1], &(xy[2]));

      integral += IntegratePolynomialsTriangle(xy, polys, order);
    }
  }

  return integral;
}


/* ******************************************************************
* Integrate polynomial over cell c.
****************************************************************** */
double NumericalIntegration::IntegratePolynomialCell(int c, const Polynomial& poly)
{
  // calculate integrals of monomials centered at cell centroid
  int order = poly.order();
  Polynomial integrals(d_, order);

  for (int k = 0; k <= order; ++k) {
    IntegrateMonomialsCell(c, k, integrals);
  }

  // dot product of coefficients of two polynomials.
  Polynomial tmp = poly;
  tmp.ChangeOrigin(mesh_->cell_centroid(c));

  double value(0.0);
  for (int k = 0; k <= order; ++k) {
    double scale = MonomialNaturalScales(c, k);
    int mk = integrals.MonomialSet(k).NumRows();
    for (int i = 0; i < mk; ++i) value += integrals(k, i) * tmp(k, i) / scale;
  }

  return value;
}


/* ******************************************************************
* Integrate over face f the product of polynomials that may have
* different origins. 
****************************************************************** */
double NumericalIntegration::IntegratePolynomialsFace(
    int f, const std::vector<const Polynomial*>& polys) const
{
  AmanziGeometry::Point enormal(d_), x1(d_), x2(d_);

  if (d_ == 2) {
    Entity_ID_List nodes;
    mesh_->face_get_nodes(f, &nodes);

    mesh_->node_get_coordinates(nodes[0], &x1);
    mesh_->node_get_coordinates(nodes[1], &x2);
    return IntegratePolynomialsEdge(x1, x2, polys);
  }

  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point fnormal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  fnormal /= area;

  // create a single polynomial centered at face centroid
  Polynomial product(d_, 0);
  product(0, 0) = 1.0;
  product.set_origin(xf);

  for (int i = 0; i < polys.size(); ++ i) {
    Polynomial tmp(*polys[i]);
    tmp.ChangeOrigin(xf);
    product *= tmp;
  }
  
  // Apply Euler theorem to each monomial
  double sum(0.0);
  Entity_ID_List edges;
  std::vector<int> dirs;

  mesh_->face_get_edges_and_dirs(f, &edges, &dirs);
  int nedges = edges.size();

  for (int n = 0; n < nedges; ++n) {
    int e = edges[n];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double length = mesh_->edge_length(e);

    enormal = tau^fnormal;

    // rescale polynomial coefficients
    double tmp = dirs[n] * ((xe - xf) * enormal) / length;

    Polynomial q(product);
    for (auto it = q.begin(); it.end() <= q.end(); ++it) {
      int m = it.MonomialSetOrder();
      int k = it.MonomialSetPosition();
      q(m, k) *= tmp / (m + 2);
    }

    // integrate along edge
    int n0, n1; 
    mesh_->edge_get_nodes(e, &n0, &n1);
    mesh_->node_get_coordinates(n0, &x1);
    mesh_->node_get_coordinates(n1, &x2);

    std::vector<const Polynomial*> q_ptr(1, &q);
    sum += IntegratePolynomialsEdge(x1, x2, q_ptr);
  }

  return sum;
}


/* ******************************************************************
* Integrate over edge (x1,x2) a product of polynomials that may have
* different origins.
****************************************************************** */
double NumericalIntegration::IntegratePolynomialsEdge(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<const Polynomial*>& polys) const
{
  // minimal quadrature rule
  int k(0);
  for (int i = 0; i < polys.size(); ++i) {
    k += polys[i]->order();
  }
  int m = k / 2;
  ASSERT(m < 8);

  AmanziGeometry::Point xm(d_);

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
* Integrate a product of polynomials that may have different origins
* over a triangle.
****************************************************************** */
double NumericalIntegration::IntegratePolynomialsTriangle(
    const std::vector<AmanziGeometry::Point>& xy,
    const std::vector<const Polynomial*>& polys, int order) const
{
  // calculate minimal quadrature rule 
  int m(order);
  if (m < 0) { 
    m = 0;
    for (int i = 0; i < polys.size(); ++i) {
      m += polys[i]->order();
    }
  }
  ASSERT(m < 10);

  int n1 = q2d_order[m][1];
  int n2 = n1 + q2d_order[m][0];

  AmanziGeometry::Point ym(d_);
  AmanziGeometry::Point y1 = xy[1] - xy[0];
  AmanziGeometry::Point y2 = xy[2] - xy[0];

  double integral(0.0);
  for (int n = n1; n < n2; ++n) { 
    ym = xy[0] + y1 * q2d_points[n][1] + y2 * q2d_points[n][2];

    double a(q2d_weights[n]);
    for (int i = 0; i < polys.size(); ++i) {
      a *= polys[i]->Value(ym);
    }
    integral += a;      
  }

  ym = y1^y2;
  double area = std::fabs(ym[0] / 2);

  return integral * area;
}


/* ******************************************************************
* Integrate over cell c a group of uniformly normalized monomials of
* the same order centered at the centroid of c.
****************************************************************** */
void NumericalIntegration::UpdateMonomialIntegralsCell(
    int c, int order, PolynomialOnMesh& integrals)
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

    for (int k = k0 + 1; k <= order; ++k) {
      IntegrateMonomialsCell(c, k, poly);
    }
  }
}


/* ******************************************************************
* Integrate over cell c a group of uniformly normalized monomials of
* the same order centered at the centroid of c: 
*    m(x) = (x - x0)^k / h^k,
* where h is a measure of cell size.
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsCell(int c, int k, Polynomial& integrals)
{
  for (int i = 0; i < integrals.MonomialSet(k).NumRows(); ++i) {
    integrals(k, i) = 0.0;
  }

  Entity_ID_List faces, nodes;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double factor = MonomialNaturalScales(c, k);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double tmp = factor * dirs[n] * ((xf - xc) * normal) / (k + d_);
    
    if (d_ == 3) {
      tmp /= mesh_->face_area(f);
      IntegrateMonomialsFace_(c, f, tmp, k, integrals);
    } else if (d_ == 2) {
      mesh_->face_get_nodes(f, &nodes);

      AmanziGeometry::Point x1(d_), x2(d_);
      mesh_->node_get_coordinates(nodes[0], &x1);
      mesh_->node_get_coordinates(nodes[1], &x2);

      x1 -= xc;  // simple change of origin
      x2 -= xc;
      IntegrateMonomialsEdge_(x1, x2, tmp, k, integrals);
    }
  }
}

/* ******************************************************************
* Integrate over face f a group of uniformly normalized monomials of
* the same order k centered at the centroid of cell c.
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsFace_(
    int c, int f, double factor, int k, Polynomial& integrals)
{
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point normal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  normal /= area;

  AmanziGeometry::Point fnormal(d_), x1(d_), x2(d_);
  std::vector<const Polynomial*> polys(1);

  PolynomialIterator it(d_);
  for (it.begin(k); it.end() <= k; ++it) {
    int l = it.MonomialSetPosition();

    // using monomial centered at xc, create polynomial centred at xf
    const int* idx = it.multi_index();
    // Polynomial poly(d_, idx, 1.0);
    // poly.set_origin(xc);
    // poly.ChangeOrigin(xf);

    Monomial mono(d_, idx, 1.0);
    mono.set_origin(xc);
    Polynomial poly = integrals.ChangeOrigin(mono, xf);

    Entity_ID_List edges;
    std::vector<int> dirs;

    mesh_->face_get_edges_and_dirs(f, &edges, &dirs);
    int nedges = edges.size();

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
      const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
      double length = mesh_->edge_length(e);

      fnormal = tau^normal;

      // rescale polynomial coefficients
      double tmp = (factor * dirs[n]) * ((xe - xf) * fnormal) / length;

      Polynomial q(poly);
      for (auto jt = poly.begin(); jt.end() <= poly.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int k = jt.MonomialSetPosition();
        q(m, k) *= tmp / (m + d_ - 1);
      }

      // integrate along edge
      int n0, n1; 
      mesh_->edge_get_nodes(e, &n0, &n1);
      mesh_->node_get_coordinates(n0, &x1);
      mesh_->node_get_coordinates(n1, &x2);

      polys[0] = &q;
      integrals(k, l) += IntegratePolynomialsEdge(x1, x2, polys);
    }
  }
}


/* ******************************************************************
* Integrate over edge (x1,x2) a group of non-normalized monomials of
* the same order k centered at zero. 
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    double factor, int k, Polynomial& integrals)
{
  AmanziGeometry::Point xm(d_);

  // minimal quadrature rule
  int m = k / 2;
  ASSERT(m < 8);

  PolynomialIterator it(d_);
  for (it.begin(k); it.end() <= k; ++it) {
    const int* idx = it.multi_index();
    int l = it.MonomialSetPosition();

    for (int n = 0; n <= m; ++n) { 
      xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);

      double a1(factor);
      for (int i = 0; i < d_; ++i) {
        a1 *= std::pow(xm[i], idx[i]);
      }

      integrals(k, l) += a1 * q1d_weights[m][n];      
    }
  }
}


/* ******************************************************************
* Pseudo-maximum value of a polynomial.
****************************************************************** */
double NumericalIntegration::PolynomialMaxValue(int f, const Polynomial& poly)
{
  int k = poly.order();
  int m = k / 2;
  double pmax;
  AmanziGeometry::Point x1(d_), x2(d_), xm(d_);

  Entity_ID_List nodes;
  mesh_->face_get_nodes(f, &nodes);

  if (d_ == 2) {
    mesh_->node_get_coordinates(nodes[0], &x1);
    mesh_->node_get_coordinates(nodes[1], &x2);

    pmax = std::max(fabs(poly.Value(x1)), fabs(poly.Value(x2)));
    for (int n = 0; n <= m; ++n) { 
      xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);
      pmax = std::max(pmax, fabs(poly.Value(xm)));
    }
  } else {
    pmax = -1.0e+99;
    for (int i = 0; i < nodes.size(); ++i) {
      mesh_->node_get_coordinates(nodes[i], &xm);
      pmax = std::max(pmax, fabs(poly.Value(xm)));
    }
  }

  return pmax;
}


/* ******************************************************************
* Re-scale polynomial coefficients: lazy implementation.
* Scaling factor is constant for monomials of the same order.
****************************************************************** */
double NumericalIntegration::MonomialNaturalScales(int c, int k) {
  if (! single_cell_) {
    basis_[c].Init(mesh_, c, k);
    return basis_[c].monomial_scales()[k];
  }

  single_cell_basis_.Init(mesh_, c, k);
  return single_cell_basis_.monomial_scales()[k];
}


void NumericalIntegration::ChangeBasisRegularToNatural(int c, Polynomial& p)
{
  for (int k = 0; k <= p.order(); ++k) {
    auto& mono = p.MonomialSet(k);
    mono /= MonomialNaturalScales(c, k);
  }
}


void NumericalIntegration::ChangeBasisNaturalToRegular(int c, Polynomial& p)
{
  for (int k = 0; k <= p.order(); ++k) {
    auto& mono = p.MonomialSet(k);
    mono *= MonomialNaturalScales(c, k);
  }
}


double NumericalIntegration::MonomialNaturalSingleScale_(int k, double volume) const {
  // return 1.0;
  return std::pow(volume, -(double)k / d_);
}

}  // namespace WhetStone
}  // namespace Amanzi

