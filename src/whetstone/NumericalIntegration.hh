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

#include "Point.hh"

#include "Monomial.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Quadrature1D.hh"
#include "Quadrature2D.hh"
#include "Quadrature3D.hh"
#include "SurfaceCoordinateSystem.hh"
#include "SurfaceMiniMesh.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

// non-member functions
double IntegrateFunctionsTriangle(
    const std::vector<AmanziGeometry::Point>& xy,
    const std::vector<const WhetStoneFunction*>& funcs, int order);

double IntegrateFunctionsTetrahedron(
    const std::vector<AmanziGeometry::Point>& xy,
    const std::vector<const WhetStoneFunction*>& funcs, int order);

double IntegrateFunctionsEdge(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<const WhetStoneFunction*>& funcs, int order);

double IntegratePolynomialsEdge(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<const PolynomialBase*>& polys);


template <class Mesh>
class NumericalIntegration { 
 public:
  NumericalIntegration(Teuchos::RCP<const Mesh> mesh);
  ~NumericalIntegration() {};

  // main methods
  // integrate product of functions with quadrature order 
  double IntegrateFunctionsTriangulatedCell(
      int c, const std::vector<const WhetStoneFunction*>& funcs, int order) const;

  double IntegrateFunctionsTriangulatedFace(
      int f, const std::vector<const WhetStoneFunction*>& funcs, int order) const;

  double IntegrateFunctionsEdge(
      int e, const std::vector<const WhetStoneFunction*>& funcs, int order) const {
    int v1, v2;
    AmanziGeometry::Point x1(d_), x2(d_);
    mesh_->edge_get_nodes(e, &v1, &v2);
    mesh_->node_get_coordinates(v1, &x1);
    mesh_->node_get_coordinates(v2, &x2);
    return WhetStone::IntegrateFunctionsEdge(x1, x2, funcs, order);
  }

  // integrate product of polynomials and monomials with different origins
  double IntegratePolynomialsCell(
      int c, const std::vector<const PolynomialBase*>& polys);

  double IntegratePolynomialsCell(
      int c, const std::vector<const PolynomialBase*>& polys,
      PolynomialOnMesh& integrals) const;

  double IntegratePolynomialsFace(
      int f, const std::vector<const PolynomialBase*>& polys) const;

  double IntegratePolynomialsEdge(
      int e, const std::vector<const PolynomialBase*>& polys) const {
    int v1, v2;
    AmanziGeometry::Point x1(d_), x2(d_);
    mesh_->edge_get_nodes(e, &v1, &v2);
    mesh_->node_get_coordinates(v1, &x1);
    mesh_->node_get_coordinates(v2, &x2);
    return WhetStone::IntegratePolynomialsEdge(x1, x2, polys);
  }

  // integral over a simplex 
  double IntegrateFunctionsSimplex(
      const std::vector<AmanziGeometry::Point>& xy,
      const std::vector<const WhetStoneFunction*>& funcs, int order) const {
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

  double IntegratePolynomialFace(int f, const Polynomial& poly) const {
    const std::vector<const PolynomialBase*> polys(1, &poly);
    return IntegratePolynomialsFace(f, polys);
  }

  double IntegratePolynomialEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const Polynomial& poly) const {
    const std::vector<const PolynomialBase*> polys(1, &poly);
    return WhetStone::IntegratePolynomialsEdge(x1, x2, polys);
  }

  // useful functions
  // -- integrate single function
  double IntegrateFunctionEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const WhetStoneFunction& func, int order) const {
    const std::vector<const WhetStoneFunction*> funcs(1, &func);
    return IntegrateFunctionsEdge(x1, x2, funcs, order);
  }

  // -- compute moments
  void CalculateFunctionMomentsEdge(
       int e, const WhetStoneFunction* func, int order,
       std::vector<double>& moments, int nquad) const;

  void CalculatePolynomialMomentsEdge(
       int e, const Polynomial& poly, int order,
       std::vector<double>& moments) const;

  void CalculatePolynomialMomentsFace(
       int f, const Polynomial& poly, int order,
       std::vector<double>& moments) const;

  // various bounds
  double PolynomialMaxValue(int f, const Polynomial& poly);

  // access
  int dimension() { return d_; }

 private:
  void IntegrateMonomialsFace_(
      int c, int f, double factor, int k, Polynomial& integrals) const;

  void IntegrateMonomialsFaceReduction_(
      int c, int f, double factor, int k, Polynomial& integrals) const;

  void IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      double factor, int k, Polynomial& integrals) const;

 private:
  Teuchos::RCP<const Mesh> mesh_;
  int d_;

  // often, we need to perform multiple integrations on the same object.
  // saving database of monomial integrals reduces computational cost, 
  // especially in 3D
  PolynomialOnMesh integrals_; 
};


//  -- trace of polynomials on manifold
Polynomial ConvertPolynomialsToSurfacePolynomial(
    const AmanziGeometry::Point& xf, 
    const std::shared_ptr<SurfaceCoordinateSystem>& coordsys,
    const std::vector<const PolynomialBase*>& polys);


/* ******************************************************************
* Constructor.
****************************************************************** */
template <class Mesh>
NumericalIntegration<Mesh>::NumericalIntegration(Teuchos::RCP<const Mesh> mesh)
  : mesh_(mesh),
    d_(mesh->space_dimension())
{};


/* ******************************************************************
* Integrate over triangulated cell a product of functions.
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::IntegrateFunctionsTriangulatedCell(
    int c, const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  double integral(0.0);

  AmanziMesh::Entity_ID_List faces, nodes;
  std::vector<AmanziGeometry::Point> xy(d_ + 1); 

  mesh_->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  xy[0] = mesh_->cell_centroid(c);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    mesh_->face_get_nodes(f, &nodes);
    int nnodes = nodes.size();

    if (d_ == 3) {
      xy[1] = mesh_->face_centroid(f);

      for (int k = 0; k < nnodes; ++k) {
        int l = (k + 1) % nnodes;
        mesh_->node_get_coordinates(nodes[k], &(xy[2]));
        mesh_->node_get_coordinates(nodes[l], &(xy[3]));

        integral += IntegrateFunctionsTetrahedron(xy, funcs, order);
      }
    } else if (d_ == 2) {
      mesh_->node_get_coordinates(nodes[0], &(xy[1]));
      mesh_->node_get_coordinates(nodes[1], &(xy[2]));

      integral += IntegrateFunctionsTriangle(xy, funcs, order);
    }
  }

  return integral;
}


/* ******************************************************************
* Integrate over triangulated face a product of functions.
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::IntegrateFunctionsTriangulatedFace(
    int f, const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  double integral(0.0);
  AmanziGeometry::Point x1(d_), x2(d_);
  AmanziMesh::Entity_ID_List faces, nodes;

  if (d_ == 3) {
    std::vector<AmanziGeometry::Point> xy(d_ + 1); 

    mesh_->face_get_nodes(f, &nodes);
    int nnodes = nodes.size();

    xy[0] = mesh_->face_centroid(f);

    for (int k = 0; k < nnodes; ++k) {
      int l = (k + 1) % nnodes;
      mesh_->node_get_coordinates(nodes[k], &(xy[1]));
      mesh_->node_get_coordinates(nodes[l], &(xy[2]));

      integral += IntegrateFunctionsTriangle(xy, funcs, order);
    }
  } else if (d_ == 2) {
    mesh_->face_get_nodes(f, &nodes);

    mesh_->node_get_coordinates(nodes[0], &x1);
    mesh_->node_get_coordinates(nodes[1], &x2);

    integral += WhetStone::IntegrateFunctionsEdge(x1, x2, funcs, order);
  }

  return integral;
}


/* ******************************************************************
* Integrate polynomial over cell c.
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::IntegratePolynomialCell(int c, const Polynomial& poly)
{
  // create/update integrals of monomials centered at cell centroid
  int order = poly.order();
  UpdateMonomialIntegralsCell(c, order, integrals_);

  // dot product of coefficients of two polynomials.
  Polynomial tmp = poly;
  tmp.ChangeOrigin(mesh_->cell_centroid(c));

  double value(0.0);
  for (int n = 0; n < tmp.size(); ++n) {
    value += integrals_.poly()(n) * tmp(n);
  }

  return value;
}


/* ******************************************************************
* Integrate product of polynomials and monomials over cells c. They 
* may have different origins.
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::IntegratePolynomialsCell(
    int c, const std::vector<const PolynomialBase*>& polys)
{
  // create a single polynomial centered at cell centroid
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  Polynomial product(d_, 0);
  product(0) = 1.0;
  product.set_origin(xc);

  for (int i = 0; i < polys.size(); ++ i) {
    Polynomial tmp(d_, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->get_origin());
    tmp.ChangeOrigin(xc);
    product *= tmp;
  }
  
  // calculate integrals of monomials centered at cell centroid
  int order = product.order();
  UpdateMonomialIntegralsCell(c, order, integrals_);

  // dot product of coefficients of two polynomials.
  double value(0.0);
  for (int n = 0; n < product.size(); ++n) {
    value += integrals_.poly()(n) * product(n);
  }

  return value;
}


/* ******************************************************************
* Integrate product of polynomials and monomials over cells c. They 
* may have different origins. Database of monomial is extended
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::IntegratePolynomialsCell(
    int c, const std::vector<const PolynomialBase*>& polys,
    PolynomialOnMesh& integrals) const
{
  // create a single polynomial centered at cell centroid
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  Polynomial product(d_, 0);
  product(0) = 1.0;
  product.set_origin(xc);

  for (int i = 0; i < polys.size(); ++ i) {
    Polynomial tmp(d_, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->get_origin());
    tmp.ChangeOrigin(xc);
    product *= tmp;
  }
  
  // calculate integrals of monomials centered at cell centroid
  int order = product.order();
  UpdateMonomialIntegralsCell(c, order, integrals);

  // dot product of coefficients of two polynomials.
  double value(0.0);
  for (int n = 0; n < product.size(); ++n) {
    value += integrals.poly()(n) * product(n);
  }

  return value;
}


/* ******************************************************************
* Integrate over face f the product of polynomials and monomials that
* may have different origins. 
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::IntegratePolynomialsFace(
    int f, const std::vector<const PolynomialBase*>& polys) const
{
  AmanziGeometry::Point enormal(d_), x1(d_), x2(d_);

  if (d_ == 2) {
    Entity_ID_List nodes;
    mesh_->face_get_nodes(f, &nodes);

    mesh_->node_get_coordinates(nodes[0], &x1);
    mesh_->node_get_coordinates(nodes[1], &x2);
    return WhetStone::IntegratePolynomialsEdge(x1, x2, polys);
  }

  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point fnormal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  fnormal /= area;

  // create a single polynomial centered at face centroid
  Polynomial product(d_, 0);
  product(0) = 1.0;
  product.set_origin(xf);

  for (int i = 0; i < polys.size(); ++ i) {
    Polynomial tmp(d_, polys[i]->order(), polys[i]->ExpandCoefficients());
    tmp.set_origin(polys[i]->get_origin());
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
    sum += WhetStone::IntegratePolynomialsEdge(x1, x2, q_ptr);
  }

  return sum;
}


/* ******************************************************************
* Integrate over cell c a group of non-normalized monomials of
* the same order centered at the centroid of c.
****************************************************************** */
template <class Mesh>
void NumericalIntegration<Mesh>::UpdateMonomialIntegralsCell(
    int c, int order, Polynomial& integrals) const
{
  int k0 = integrals.order();

  if (k0 < order) {
    integrals.Reshape(d_, order);

    for (int k = k0 + 1; k <= order; ++k)
      IntegrateMonomialsCell(c, k, integrals);
  }
}


template <class Mesh>
void NumericalIntegration<Mesh>::UpdateMonomialIntegralsCell(
    int c, int order, PolynomialOnMesh& integrals) const
{
  Polynomial& poly = integrals.poly();
  int k0 = poly.order();

  // data are uniquecly defined by topological dimension of
  // a geometric entity and its mesh id. 
  if (integrals.get_kind() != (Entity_kind)WhetStone::CELL ||
      integrals.get_id() != c ||
      poly.dimension() != d_) {
    integrals.set_kind((Entity_kind)WhetStone::CELL);
    integrals.set_id(c);
    k0 = -1;
  }

  // add additional integrals of monomials
  if (k0 < order) {
    poly.Reshape(d_, order);

    for (int k = k0 + 1; k <= order; ++k)
      IntegrateMonomialsCell(c, k, poly);
  }
}


/* ******************************************************************
* Integrate over cell c a group of non-normalized monomials of
* the same order centered at the centroid of c.
****************************************************************** */
template <class Mesh>
void NumericalIntegration<Mesh>::IntegrateMonomialsCell(
    int c, int k, Polynomial& integrals) const
{
  if (k == 0) {
    integrals(0) = mesh_->cell_volume(c);
    return;
  } 

  if (k == 1) {
    for (int i = 0; i < d_; ++i) integrals(1 + i) = 0.0;
    return;
  }

  int nk = PolynomialSpaceDimension(d_, k - 1);
  int mk = MonomialSpaceDimension(d_, k);
  for (int i = 0; i < mk; ++i) {
    integrals(nk + i) = 0.0;
  }

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
      IntegrateMonomialsFaceReduction_(c, f, tmp, k, integrals);
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
* Integrate over face f a group of non-normalized monomials of
* the same order k centered at the centroid of cell c.
****************************************************************** */
template <class Mesh>
void NumericalIntegration<Mesh>::IntegrateMonomialsFace_(
    int c, int f, double factor, int k, Polynomial& integrals) const
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
      for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int s = jt.PolynomialPosition();
        q(s) *= tmp / (m + d_ - 1);
      }

      // integrate along edge (based on Euler theorem)
      int n0, n1; 
      mesh_->edge_get_nodes(e, &n0, &n1);
      mesh_->node_get_coordinates(n0, &x1);
      mesh_->node_get_coordinates(n1, &x2);

      polys[0] = &q;
      integrals(nk + l) += WhetStone::IntegratePolynomialsEdge(x1, x2, polys);

      // integrate along edge (based on change of variables)
      /*
      std::vector<AmanziGeometry::Point> tau_edge(1, tau);
      q.ChangeCoordinates(xe, tau_edge);  

      int m(1);
      double sum(0.0);
      for (int i = 0; i < q.size(); i += 2) {
        sum += q(i) / (i + 1) / m;
        m *= 4;
      }
      integrals(nk + l) += sum * length;
      */
    }
  }
}


/* ******************************************************************
* Integrate over edge (x1,x2) a group of non-normalized monomials of
* the same order k centered at zero. 
****************************************************************** */
template <class Mesh>
void NumericalIntegration<Mesh>::IntegrateMonomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    double factor, int k, Polynomial& integrals) const
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
* Calculate moments upto given order.
****************************************************************** */
template <class Mesh>
void NumericalIntegration<Mesh>::CalculateFunctionMomentsEdge(
    int e, const WhetStoneFunction* func, int order,
    std::vector<double>& moments, int nquad) const
{
  moments.clear();

  double len = mesh_->edge_length(e);

  if (order >= 0) {
    std::vector<const WhetStone::WhetStoneFunction*> funcs(1);
    funcs[0] = func;

    double tmp = IntegrateFunctionsEdge(e, funcs, nquad) / len;
    moments.push_back(tmp);
  }

  if (order > 0) {
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);

    std::vector<AmanziGeometry::Point> coordsys(1, tau);
    std::vector<const WhetStone::WhetStoneFunction*> funcs(2);
    funcs[0] = func;

    PolynomialIterator it(1);
    for (it.begin(1); it.MonomialSetOrder() <= order; ++it) {
      const int* index = it.multi_index();
      Polynomial emono(1, index, 1.0);
      emono.InverseChangeCoordinates(xe, coordsys);  
      funcs[1] = &emono;

      double tmp = IntegrateFunctionsEdge(e, funcs, nquad) / len;
      moments.push_back(tmp);
    }
  }
}


/* ******************************************************************
* Calculate moments upto given order.
****************************************************************** */
template <class Mesh>
void NumericalIntegration<Mesh>::CalculatePolynomialMomentsEdge(
    int e, const Polynomial& poly, int order, std::vector<double>& moments) const
{
  moments.clear();

  int m = poly.order();
  double len = mesh_->edge_length(e);

  if (order >= 0) {
    std::vector<const WhetStone::WhetStoneFunction*> funcs(1);
    funcs[0] = &poly;

    double tmp = IntegrateFunctionsEdge(e, funcs, m) / len;
    moments.push_back(tmp);
  }

  if (order > 0) {
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);

    Polynomial q(poly);
    std::vector<AmanziGeometry::Point> coordsys(1, tau);
    q.ChangeCoordinates(xe, coordsys);  

    std::vector<const WhetStone::WhetStoneFunction*> funcs(2);
    funcs[0] = &q;

    AmanziGeometry::Point x1(1), x2(1);
    x1[0] = -0.5;
    x2[0] =  0.5;

    PolynomialIterator it(1);
    for (it.begin(1); it.MonomialSetOrder() <= order; ++it) {
      int m2 = it.MonomialSetOrder();
      const int* index = it.multi_index();
      Polynomial emono(1, index, 1.0);
      funcs[1] = &emono;

      double tmp = WhetStone::IntegrateFunctionsEdge(x1, x2, funcs, m + m2);
      moments.push_back(tmp);
    }
  }
}


/* ******************************************************************
* Approximate maximum value of a polynomial.
****************************************************************** */
template <class Mesh>
double NumericalIntegration<Mesh>::PolynomialMaxValue(int f, const Polynomial& poly)
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

}  // namespace WhetStone
}  // namespace Amanzi

#endif
