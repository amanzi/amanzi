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
  // extend list of integrals of monomials
  UpdateIntegrals_(c, 2 * order_);
  const Polynomial& integrals = integrals_[c];
   
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

  ChangeBasis_(c, M);

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

  // extend list of integrals of monomials
  int uk(Kcopy.order());
  UpdateIntegrals_(c, 2 * order_ + uk);
  const Polynomial& integrals = integrals_[c];
   
  // sum up integrals to the mass matrix
  int multi_index[3];
  double ak, bk, al, bl;
  Polynomial p(d_, order_), q(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);
  M.PutScalar(0.0);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

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

  ChangeBasis_(c, M);

  return 0;
}


/* ******************************************************************
* Advection matrix for Taylor basis and cell-based velocity uc.
****************************************************************** */
int DG_Modal::AdvectionMatrixPoly(int c, const VectorPolynomial& u, DenseMatrix& A, bool grad_on_test)
{
  // rebase the polynomial
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  VectorPolynomial ucopy(u);
  for (int i = 0; i < d_; ++i) {
    ucopy[i].ChangeOrigin(xc);
  }

  // extend list of integrals of monomials
  int uk(ucopy[0].order());
  UpdateIntegrals_(c, 2 * order_ + std::max(0, uk - 1));
  const Polynomial& integrals = integrals_[c];
   
  // sum-up integrals to the advection matrix
  int multi_index[3];
  Polynomial p(d_, order_), q(d_, order_);
  VectorPolynomial pgrad;

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

  ChangeBasis_(c, A);

  // gradient operator is applied to solution
  if (!grad_on_test) {
    A.Transpose();
  }

  return 0;
}


/* ******************************************************************
* Flux matrix for Taylor basis and normal velocity u.n.
* Velocity is given in the face-based Taylor basis.
* If jump_on_test=true, we calculate
* 
*   \Int { (u.n) \rho^* [\psi] } dS
* 
* where star means downwind, \psi is a test function, and \rho is 
* a solution. Otherwise, we calculate 
*
*   \Int { (u.n) \psi^* [\rho] } dS
****************************************************************** */
int DG_Modal::FluxMatrixPoly(int f, const Polynomial& un, DenseMatrix& A,
                             bool jump_on_test)
{
  AmanziMesh::Entity_ID_List cells, nodes;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  Polynomial poly0(d_, order_), poly1(d_, order_);
  int size = poly0.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  // identify index of downwind cell (id)
  int dir, id(0); 
  mesh_->face_normal(f, false, cells[0], &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  double vel = un.Value(xf) * dir;
  if (vel > 0.0) {
    if (ncells == 1) return 0;
    id = 1;
  } else {
    if (ncells == 1) return 0;  // Assume that u.n=0 FIXME
    dir = -dir;
  }
  int col(id * size);
  int row(size - col);

  // Calculate integrals needed for scaling
  UpdateIntegrals_(cells[id], 2 * order_);
  const Polynomial& integrals_dw = integrals_[cells[id]];

  UpdateIntegrals_(cells[1 - id], 2 * order_);
  const Polynomial& integrals_up = integrals_[cells[1 - id]];

  // integrate traces of polynomials on face f
  std::vector<const Polynomial*> polys(3);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);

    Polynomial p0(d_, idx0), p1(d_, idx0);
    p0.set_origin(mesh_->cell_centroid(cells[id]));
    p1.set_origin(mesh_->cell_centroid(cells[1 - id]));

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);

      Polynomial q(d_, idx1);
      q.set_origin(mesh_->cell_centroid(cells[id]));

      polys[0] = &un;
      polys[1] = &p0;
      polys[2] = &q;

      // downwind-downwind integral
      double vel1 = IntegratePolynomialsFace_(f, polys);
      vel1 /= mesh_->face_area(f);
      vel1 *= dir;  

      // upwind-downwind integral
      polys[1] = &p1;

      double vel0 = IntegratePolynomialsFace_(f, polys);
      vel0 /= mesh_->face_area(f);
      vel0 *= dir;  

      if (ncells == 1) {
        A(k, l) = vel1;
      } else {
        A(row + k, col + l) = vel0;
        A(col + k, col + l) = -vel1;
      }
    }
  }

  ChangeBasis_(cells[0], cells[1], A);

  // gradient operator is applied to solution
  if (!jump_on_test) {
    A.Transpose();
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
void DG_Modal::IntegrateMonomialsFace_(
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
      monomials.coefs()[l] += IntegratePolynomialsEdge_(x1, x2, polys);
    }
  }
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
* Integrate over face f the product of polynomials that may have
* different origin. 
****************************************************************** */
double DG_Modal::IntegratePolynomialsFace_(
    int f, const std::vector<const Polynomial*>& polys) const
{
  AmanziGeometry::Point enormal(d_), x1(d_), x2(d_);

  if (d_ == 2) {
    Entity_ID_List nodes;
    mesh_->face_get_nodes(f, &nodes);

    mesh_->node_get_coordinates(nodes[0], &x1);
    mesh_->node_get_coordinates(nodes[1], &x2);
    return IntegratePolynomialsEdge_(x1, x2, polys);
  }

  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  AmanziGeometry::Point fnormal = mesh_->face_normal(f);
  double area = mesh_->face_area(f);
  fnormal /= area;

  // create a single polynomial centered at face centroid
  Polynomial product(d_, 0);
  product(0, 0) = 1.0;
  product.ChangeOrigin(xf);

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
    sum += IntegratePolynomialsEdge_(x1, x2, q_ptr);
  }

  return sum;
}


/* ******************************************************************
* Integrate over edge (x1,x2) a product of polynomials that may have
* different origins.
****************************************************************** */
double DG_Modal::IntegratePolynomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    const std::vector<const Polynomial*>& polys) const
{
  // minimal quadrature rule
  int k(0);
  for (int i = 0; i < polys.size(); ++i) {
    k += polys[i]->order();
  }
  int m = k / 2;

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
* Change basis from unscaled monomials to Taylor basis for a cell
****************************************************************** */
void DG_Modal::ChangeBasis_(int c, DenseMatrix& A)
{
  // optional update of integrals database
  UpdateIntegrals_(c, 2 * order_);

  int nrows = A.NumRows();
  double ak, bk;
  std::vector<double> a(nrows), b(nrows);

  Iterator it(d_);
  for (it.begin(); it.end() <= order_; ++it) {
    int k = it.PolynomialPosition();
    TaylorBasis(c, it, &ak, &bk);
    a[k] = ak;
    b[k] = -ak * bk;
  }

  // calculate A * R
  for (int k = 1; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(i, k) = A(i, k) * a[k] + A(i, 0) * b[k];
    }
  }

  // calculate R^T * A * R
  for (int k = 1; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(k, i) = A(k, i) * a[k] + A(0, i) * b[k];
    }
  }
}


/* ******************************************************************
* Change basis from unscaled monomials to Taylor basis for a face
****************************************************************** */
void DG_Modal::ChangeBasis_(int c1, int c2, DenseMatrix& A)
{
  // optional update of integrals database
  UpdateIntegrals_(c1, 2 * order_);
  UpdateIntegrals_(c2, 2 * order_);

  int m, nrows = A.NumRows();
  double ak, bk;
  std::vector<double> a(nrows), b(nrows);

  m = nrows / 2;
  Iterator it(d_);
  for (it.begin(); it.end() <= order_; ++it) {
    int k = it.PolynomialPosition();

    TaylorBasis(c1, it, &ak, &bk);
    a[k] = ak;
    b[k] = -ak * bk;

    TaylorBasis(c2, it, &ak, &bk);
    a[m + k] = ak;
    b[m + k] = -ak * bk;
  }

  // calculate A * R
  for (int k = 1; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(i, k) = A(i, k) * a[k] + A(i, 0) * b[k];
    }
  }

  // calculate R^T * A * R
  for (int k = 1; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(k, i) = A(k, i) * a[k] + A(0, i) * b[k];
    }
  }
}


/* ******************************************************************
* Transform monomial \psi_k to polynomial a (\psi_k - b \psi_0) where
* factor b orthogonolizes to constant and factor normalizes to 1.
* NOTE: Both polynomilas are centered at cell centroid.
****************************************************************** */
void DG_Modal::TaylorBasis(int c, const Iterator& it, double* a, double* b)
{
  const int* multi_index = it.multi_index();

  int n(0);
  for (int i = 0; i < d_; ++i) {
    n += multi_index[i];
  }

  // We do not modify the first function
  if (n == 0) {
    *a = 1.0;
    *b = 0.0;
  } else if (basis_ == TAYLOR_BASIS_SIMPLE) {
    *a = 1.0;
    *b = 0.0;
  } else {
    UpdateScales_(c, order_);
    *a = scales_a_[c].monomials(n).coefs()[it.MonomialPosition()];
    *b = scales_b_[c].monomials(n).coefs()[it.MonomialPosition()];
  }
}


/* ******************************************************************
* Calculate polynomial using basis and coefficients.
****************************************************************** */
Polynomial DG_Modal::CalculatePolynomial(int c, const std::vector<double>& coefs) const
{
  ASSERT(scales_a_.size() != 0);
  ASSERT(scales_a_[c].size() == coefs.size());

  Polynomial poly(scales_a_[c]);
  poly.set_origin(mesh_->cell_centroid(c));

  auto jt = coefs.begin();
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    poly(m, k) *= *jt; 
    poly(0, 0) -= *jt * scales_b_[c](m, k);
    ++jt;
  }

  return poly;
}


/* ******************************************************************
* Update integrals of non-normalized monomials.
****************************************************************** */
void DG_Modal::UpdateIntegrals_(int c, int order)
{
  if (integrals_.size() == 0) {
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    integrals_.resize(ncells_wghost);

    for (int n = 0; n < ncells_wghost; ++n) {
      integrals_[n].Reshape(d_, 0);
      integrals_[n](0, 0) = mesh_->cell_volume(n);
    }
  }

  // add additional integrals of monomials
  int k0 = integrals_[c].order();
  if (k0 < order) {
    integrals_[c].Reshape(d_, order);

    for (int k = k0 + 1; k <= order; ++k) {
      IntegrateMonomialsCell_(c, integrals_[c].monomials(k));
    }
  }
}


/* ******************************************************************
* Partial orthonormalization of Taylor basis functions.
****************************************************************** */
void DG_Modal::UpdateScales_(int c, int order)
{
  if (scales_a_.size() == 0) {
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    scales_a_.resize(ncells_wghost);
    scales_b_.resize(ncells_wghost);

    for (int n = 0; n < ncells_wghost; ++n) {
      scales_a_[n].Reshape(d_, order);
      scales_b_[n].Reshape(d_, order);
    }

    // For the moment, we update everything in one shot
    for (int n = 0; n < ncells_wghost; ++n) {
      UpdateIntegrals_(n, 2 * order);

      const Polynomial& integrals = integrals_[n];
      Polynomial poly(d_, order);

      double a, b, norm;
      double volume = integrals(0, 0); 

      for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
        int k = it.MonomialPosition();
        const int* multi_index = it.multi_index();
        int index[d_]; 

        int m(0);
        for (int i = 0; i < d_; ++i) {
          m += multi_index[i];
          index[i] = 2 * multi_index[i];
        }

        if (m == 0) {
          a = 1.0;
          b = 0.0;
        } else {
          const auto& aux1 = integrals.monomials(m).coefs();
          b = aux1[k] / volume;

          const auto& aux2 = integrals.monomials(2 * m).coefs();
          norm = aux2[integrals.MonomialPosition(index)];
          norm -= b * b * volume;

          a = std::pow(volume / norm, 0.5);
        }

        scales_a_[n].monomials(m).coefs()[k] = a;
        scales_b_[n].monomials(m).coefs()[k] = b;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

