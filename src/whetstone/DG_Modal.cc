/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin methods.
*/

#include "Point.hh"

#include "DenseMatrix.hh"
#include "DG_Modal.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Mass matrix for Taylor basis functions. 
****************************************************************** */
int DG_Modal::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  // calculate all monomials
  Polynomial integrals(d_, 2 * order_);

  integrals.monomials(0).coefs()[0] = mesh_->cell_volume(c);

  for (int k = 1; k <= integrals.order(); ++k) {
    IntegrateMonomialsCell_(c, integrals.monomials(k));
  }
   
  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_), q(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = p.PolynomialPosition(idx_p);

    for (auto jt = q.begin(); jt.end() <= q.end(); ++jt) {
      const int* idx_q = jt.multi_index();
      int l = q.PolynomialPosition(idx_q);
      
      int n(0);
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = idx_p[i] + idx_q[i];
        n += multi_index[i];
      }

      const auto& coefs = integrals.monomials(n).coefs();
      M(k, l) = K(0, 0) * coefs[p.MonomialPosition(multi_index)];
    }
  }

  return 0;
}


/* ******************************************************************
* Mass matrix for Taylor basis and polynomial coefficient K.
****************************************************************** */
int DG_Modal::MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M)
{
  // rebase the polynomial
  Polynomial Kcopy(K);
  Kcopy.ChangeOrigin(mesh_->cell_centroid(c));

  // calculate integrals of monomials centered at cell centroid
  int uk(Kcopy.order());
  Polynomial integrals(d_, 2 * order_ + uk);

  integrals.monomials(0).coefs()[0] = mesh_->cell_volume(c);

  for (int k = 1; k <= integrals.order(); ++k) {
    IntegrateMonomialsCell_(c, integrals.monomials(k));
  }

  // sum up integrals to the mass matrix
  int multi_index[3];
  Polynomial p(d_, order_), q(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);
  M.PutScalar(0.0);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = p.PolynomialPosition(idx_p);

    for (auto mt = Kcopy.begin(); mt.end() <= Kcopy.end(); ++mt) {
      const int* idx_K = mt.multi_index();
      int m = mt.MonomialPosition();
      double factor = Kcopy.monomials(mt.end()).coefs()[m];
      if (factor == 0.0) continue;

      for (auto jt = q.begin(); jt.end() <= q.end(); ++jt) {
        const int* idx_q = jt.multi_index();
        int l = q.PolynomialPosition(idx_q);

        int n(0);
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = idx_p[i] + idx_q[i] + idx_K[i];
          n += multi_index[i];
        }

        const auto& coefs = integrals.monomials(n).coefs();
        M(k, l) += factor * coefs[p.MonomialPosition(multi_index)];
      }
    }
  }

  return 0;
}


/* ******************************************************************
* Advection matrix for Taylor basis and cell-based velocity uc.
****************************************************************** */
int DG_Modal::AdvectionMatrixPoly(int c, const VectorPolynomial& u, DenseMatrix& A)
{
  // rebase the polynomial
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  VectorPolynomial ucopy(u);
  for (int i = 0; i < d_; ++i) {
    ucopy[i].ChangeOrigin(xc);
  }

  // calculate integrals of monomials centered at cell centroid
  int uk(ucopy[0].order());
  Polynomial integrals(d_, 2 * order_ + uk - 1);
  VectorPolynomial pgrad;

  integrals.monomials(0).coefs()[0] = 0.0;

  for (int k = 1; k <= integrals.order(); ++k) {
    IntegrateMonomialsCell_(c, integrals.monomials(k));
  }

  // sum-up integrals to the advection matrix
  int multi_index[3];
  Polynomial p(d_, order_), q(d_, order_);

  int nrows = p.size();
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = p.PolynomialPosition(idx_p);

    // product of polynomials need to align origins
    Polynomial pp(d_, idx_p);
    pp.set_origin(xc);

    pp.Gradient(pgrad);
    Polynomial tmp(pgrad * ucopy);

    for (auto mt = tmp.begin(); mt.end() <= tmp.end(); ++mt) {
      const int* idx_K = mt.multi_index();
      int m = mt.MonomialPosition();
      double factor = tmp.monomials(mt.end()).coefs()[m];
      if (factor == 0.0) continue;

      for (auto jt = q.begin(); jt.end() <= q.end(); ++jt) {
        const int* idx_q = jt.multi_index();
        int l = q.PolynomialPosition(idx_q);

        int n(0);
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = idx_q[i] + idx_K[i];
          n += multi_index[i];
        }

        const auto& coefs = integrals.monomials(n).coefs();
        A(k, l) += factor * coefs[p.MonomialPosition(multi_index)];
      }
    }
  }

  if (nrows > 1) std::cout << A << std::endl;

  return 0;
}


/* ******************************************************************
* Flux matrix for Taylor basis and normal velocity u.n.
* Velocity is given in the face-based Taylor basis.
****************************************************************** */
int DG_Modal::FluxMatrixPoly(int f, const Polynomial& un, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells, nodes;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  Polynomial poly0(d_, order_), poly1(d_, order_);
  int size = poly0.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  // identify downwind cell
  int dir, id(0); 
  mesh_->face_normal(f, false, cells[0], &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  double vel = un.Value(xf) * dir;
  if (vel > 0.0) {
    if (ncells == 1) return 0;
    id = 1;
  } else {
    dir = -dir;
  }
  int col(id * size);
  int row(size - col);

  // integrate traces of polynomials on face f
  mesh_->face_get_nodes(f, &nodes);

  AmanziGeometry::Point x1(d_), x2(d_);
  mesh_->node_get_coordinates(nodes[0], &x1);
  mesh_->node_get_coordinates(nodes[1], &x2);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);

    Polynomial p(d_, idx0);
    p.set_origin(mesh_->cell_centroid(cells[id]));

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);

      Polynomial q(d_, idx1);
      q.set_origin(mesh_->cell_centroid(cells[id]));

      std::vector<Polynomial> polys;
      polys.push_back(un);
      polys.push_back(p);
      polys.push_back(q);

      vel = IntegratePolynomialsEdge_(x1, x2, polys);
      vel /= mesh_->face_area(f);
      vel *= dir;  

      if (ncells == 1) {
        A(k, l) = vel;
      } else {
        A(row + k, col + l) = vel;
        A(col + k, col + l) = -vel;
      }
    }
  }

  return 0;
}


/* ******************************************************************
* Integrate over cell c a group of non-normalized monomials of the
* same order centered at the centroid of c.
****************************************************************** */
void DG_Modal::IntegrateMonomialsCell_(int c, Monomial& monomials)
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
      IntegrateMonomialsFace_(f, tmp, monomials);
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
void DG_Modal::IntegrateMonomialsFace_(int f, double factor, Monomial& monomials)
{
}


/* ******************************************************************
* Integrate over edge (x1,x2) a group of non-normalized monomials of
* the same order centered at zero. 
****************************************************************** */
void DG_Modal::IntegrateMonomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    double factor, Monomial& monomials)
{
  AmanziGeometry::Point xm(d_);
  auto& coefs = monomials.coefs();

  // minimal quadrature rule
  int k = monomials.order();
  int m = k / 2;

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


/* ******************************************************************
* Integrate over edge (x1,x2) a product of polynomials that may have
* different origins.
****************************************************************** */
double DG_Modal::IntegratePolynomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<Polynomial>& polys) const
{
  // minimal quadrature rule
  int k(0);
  for (int i = 0; i < polys.size(); ++i) {
    k += polys[i].order();
  }
  int m = k / 2;

  AmanziGeometry::Point xm(d_);

  double integral(0.0);
  for (int n = 0; n <= m; ++n) { 
    xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);

    double a(q1d_weights[m][n]);
    for (int i = 0; i < polys.size(); ++i) {
      a *= polys[i].Value(xm);
    }
    integral += a;      
  }

  return integral * norm(x2 - x1);
}

}  // namespace WhetStone
}  // namespace Amanzi

