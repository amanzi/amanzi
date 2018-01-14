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
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Integrate polynomial over cell c.
****************************************************************** */
double NumericalIntegration::IntegratePolynomialCell(int c, const Polynomial& poly)
{
  // calculate integrals of monomials centered at cell centroid
  int order = poly.order();
  Polynomial integrals(d_, order);

  for (int k = 0; k <= order; ++k) {
    IntegrateMonomialsCell(c, integrals.monomials(k));
  }

  // dot product of coefficients of two polynomials.
  Polynomial tmp = poly;
  tmp.ChangeOrigin(mesh_->cell_centroid(c));

  double value(0.0);
  for (int k = 0; k <= order; ++k) {
    std::vector<double>& val1 = integrals.monomials(k).coefs();
    std::vector<double>& val2 = tmp.monomials(k).coefs();
    for (int i = 0; i < val1.size(); ++i) value += val1[i] * val2[i];
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
      int m = it.MonomialOrder();
      int k = it.MonomialPosition();
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
  ASSERT(m < 5);

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
* Integrate over cell c a group of non-normalized monomials of the
* same order centered at the centroid of c.
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsCell(int c, Monomial& monomials)
{
  int k = monomials.order();

  Entity_ID_List faces, nodes;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double tmp = dirs[n] * ((xf - xc) * normal) / (k + d_);
    
    if (d_ == 3) {
      tmp /= mesh_->face_area(f);
      IntegrateMonomialsFace_(c, f, tmp, monomials);
    } else if (d_ == 2) {
      mesh_->face_get_nodes(f, &nodes);

      AmanziGeometry::Point x1(d_), x2(d_);
      mesh_->node_get_coordinates(nodes[0], &x1);
      mesh_->node_get_coordinates(nodes[1], &x2);

      x1 -= xc;  // simple change of origin
      x2 -= xc;
      IntegrateMonomialsEdge_(x1, x2, tmp, monomials);
    }
  }
}


/* ******************************************************************
* Integrate over face f a group of non-normalized monomials of the
* same order centered at the centroid of cell c.
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsFace_(
    int c, int f, double factor, Monomial& monomials)
{
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point normal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  normal /= area;

  AmanziGeometry::Point fnormal(d_), x1(d_), x2(d_);
  std::vector<const Polynomial*> polys(1);

  int k = monomials.order();
  for (auto it = monomials.begin(); it.end() <= k; ++it) {
    int l = it.MonomialPosition();

    // using monomial centered at xc, create polynomial centred at xf
    const int* idx = it.multi_index();
    Polynomial poly(d_, idx);
    poly.set_origin(xc);
    poly.ChangeOrigin(xf);

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
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        q(m, k) *= tmp / (m + d_ - 1);
      }

      // integrate elong edge
      int n0, n1; 
      mesh_->edge_get_nodes(e, &n0, &n1);
      mesh_->node_get_coordinates(n0, &x1);
      mesh_->node_get_coordinates(n1, &x2);

      polys[0] = &q;
      monomials.coefs()[l] += IntegratePolynomialsEdge(x1, x2, polys);
    }
  }
}


/* ******************************************************************
* Integrate over edge (x1,x2) a group of non-normalized monomials of
* the same order centered at zero. 
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    double factor, Monomial& monomials)
{
  AmanziGeometry::Point xm(d_);
  auto& coefs = monomials.coefs();

  // minimal quadrature rule
  int k = monomials.order();
  int m = k / 2;
  ASSERT(m < 5);

  for (auto it = monomials.begin(); it.end() <= k; ++it) {
    const int* idx = it.multi_index();
    int l = it.MonomialPosition();

    for (int n = 0; n <= m; ++n) { 
      xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);

      double a1(factor);
      for (int i = 0; i < d_; ++i) {
        a1 *= std::pow(xm[i], idx[i]);
      }

      coefs[l] += a1 * q1d_weights[m][n];      
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

