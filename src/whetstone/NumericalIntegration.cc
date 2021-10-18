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
#include "SingleFaceMesh.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor.
****************************************************************** */
NumericalIntegration::NumericalIntegration(Teuchos::RCP<const AmanziMesh::MeshLight> mesh)
  : mesh_(mesh),
    d_(mesh->space_dimension())
{};


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
      a *= funcs[i]->Value(ym);
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
      a *= funcs[i]->Value(ym);
    }
    integral += a;      
  }

  double volume = std::fabs(((y1^y2) * y3) / 6);

  return integral * volume;
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
      a *= funcs[i]->Value(xm);
    }
    integral += a;      
  }

  return integral * norm(x2 - x1);
}


/* ******************************************************************
* Integrate over triangulated cell a product of functions.
****************************************************************** */
double NumericalIntegration::IntegrateFunctionsTriangulatedCell(
    int c, const std::vector<const WhetStoneFunction*>& funcs, int order) const
{
  double integral(0.0);

  AmanziMesh::Entity_ID_List nodes;
  std::vector<AmanziGeometry::Point> xy(d_ + 1); 

  const auto& faces = mesh_->cell_get_faces(c);
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
double NumericalIntegration::IntegrateFunctionsTriangulatedFace(
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
double NumericalIntegration::IntegratePolynomialCell(int c, const Polynomial& poly)
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
double NumericalIntegration::IntegratePolynomialsCell(
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
double NumericalIntegration::IntegratePolynomialsCell(
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
double NumericalIntegration::IntegratePolynomialsFace(
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
* Integrate over cell c a group of non-normalized monomials of
* the same order centered at the centroid of c.
****************************************************************** */
void NumericalIntegration::UpdateMonomialIntegralsCell(
    int c, int order, Polynomial& integrals) const
{
  int k0 = integrals.order();

  if (k0 < order) {
    integrals.Reshape(d_, order);

    for (int k = k0 + 1; k <= order; ++k)
      IntegrateMonomialsCell(c, k, integrals);
  }
}


void NumericalIntegration::UpdateMonomialIntegralsCell(
    int c, int order, PolynomialOnMesh& integrals) const
{
  Polynomial& poly = integrals.poly();
  int k0 = poly.order();

  // reset polynomial metadata
  if (integrals.get_kind() != (Entity_kind)WhetStone::CELL || integrals.get_id() != c) {
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
void NumericalIntegration::IntegrateMonomialsCell(
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

  Entity_ID_List nodes;

  // mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  const auto& faces = mesh_->cell_get_faces(c);
  const auto& dirs = mesh_->cell_get_face_dirs(c);
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
void NumericalIntegration::IntegrateMonomialsFace_(
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
void NumericalIntegration::IntegrateMonomialsEdge_(
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
* Approximate maximum value of a polynomial.
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


/* ******************************************************************
* Integrate over face f a group of non-normalized monomials of
* the same order k centered at the centroid of cell c.
****************************************************************** */
void NumericalIntegration::IntegrateMonomialsFaceReduction_(
    int c, int f, double factor, int k, Polynomial& integrals) const
{
  int nk = PolynomialSpaceDimension(d_, k - 1);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);

  // create a surface mesh
  SurfaceCoordinateSystem coordsys(xf, normal);
  Teuchos::RCP<const SingleFaceMesh> surf_mesh = Teuchos::rcp(new SingleFaceMesh(mesh_, f, coordsys));
  NumericalIntegration numi_f(surf_mesh);

  PolynomialIterator it(d_);
  for (it.begin(k); it.MonomialSetOrder() <= k; ++it) {
    int l = it.MonomialSetPosition();

    // using monomial centered at xc, create 2D polynomial centered at xf
    const int* idx = it.multi_index();
    Polynomial poly(d_, idx, 1.0);
    poly.set_origin(xc);
    poly.ChangeCoordinates(xf, *coordsys.tau());  

    integrals(nk + l) += factor * numi_f.IntegratePolynomialCell(0, poly);
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

