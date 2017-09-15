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
    DenseMatrix Atmp(A);
    A.Transpose(Atmp);
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
  mesh_->face_get_nodes(f, &nodes);

  AmanziGeometry::Point x1(d_), x2(d_);
  mesh_->node_get_coordinates(nodes[0], &x1);
  mesh_->node_get_coordinates(nodes[1], &x2);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);

    Polynomial p(d_, idx0);
    p.set_origin(mesh_->cell_centroid(cells[id]));

    // scaling of basis functions
    double ak0, bk0, ak1, bk1, al1, bl1;
    TaylorBasis(cells[id], it, &ak1, &bk1);
    TaylorBasis(cells[1 - id], it, &ak0, &bk0);

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);

      Polynomial q(d_, idx1);
      q.set_origin(mesh_->cell_centroid(cells[id]));

      std::vector<Polynomial> polys;
      polys.push_back(un);
      polys.push_back(p);
      polys.push_back(q);

      // downwind-downwind integral
      double vel1 = IntegratePolynomialsEdge_(x1, x2, polys);
      vel1 /= mesh_->face_area(f);
      vel1 *= dir;  

      TaylorBasis(cells[id], jt, &al1, &bl1);
      vel1 *= ak1 * al1;  // FIXME (up to linear basis functions)

      // upwind-downwind integral
      polys[1].set_origin(mesh_->cell_centroid(cells[1 - id]));

      double vel0 = IntegratePolynomialsEdge_(x1, x2, polys);
      vel0 /= mesh_->face_area(f);
      vel0 *= dir;  
      vel0 *= ak0 * al1;

      if (ncells == 1) {
        A(k, l) = vel1;
      } else {
        A(row + k, col + l) = vel0;
        A(col + k, col + l) = -vel1;
      }
    }
  }

  // gradient operator is applied to solution
  if (!jump_on_test) {
    DenseMatrix Atmp(A);
    A.Transpose(Atmp);
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


/* ******************************************************************
* Change basis from unscaled monomials to Taylor basis
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
* Update integrals of non-normalized monomials.
****************************************************************** */
void DG_Modal::UpdateIntegrals_(int c, int order)
{
  if (integrals_.size() == 0) {
    int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    integrals_.resize(ncells_owned);

    for (int n = 0; n < ncells_owned; ++n) {
      integrals_[n].Reshape(d_, 0);
      integrals_[n].monomials(0).coefs()[0] = mesh_->cell_volume(n);
    }
  }

  // add optionally additional integrals of monomials
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
    int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    scales_a_.resize(ncells_owned);
    scales_b_.resize(ncells_owned);

    for (int n = 0; n < ncells_owned; ++n) {
      scales_a_[n].Reshape(d_, order);
      scales_b_[n].Reshape(d_, order);
    }

    // For the moment, we update everything in one shot
    for (int n = 0; n < ncells_owned; ++n) {
      UpdateIntegrals_(n, 2 * order);

      const Polynomial& integrals = integrals_[n];
      Polynomial poly(d_, order);

      double a, b, norm;
      double volume = integrals.monomials(0).coefs()[0]; 

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

          scales_a_[n].monomials(m).coefs()[k] = a;
          scales_b_[n].monomials(m).coefs()[k] = b;
        }
      }
    }
  }
}


/* ******************************************************************
* Error calculation requires geometric center.
****************************************************************** */
AmanziGeometry::Point DG_Modal::cell_geometric_center(int c) const
{
  Entity_ID_List nodes;
  AmanziGeometry::Point v(d_), xg(d_);

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();
  for (int i = 0; i < nnodes; ++i) {
    mesh_->node_get_coordinates(nodes[i], &v);
    xg += v;
  } 
  xg /= nnodes;
  
  return xg;
}

}  // namespace WhetStone
}  // namespace Amanzi

