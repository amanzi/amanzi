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
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Mass matrix for Taylor basis functions. 
****************************************************************** */
int DG_Modal::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  double K00 = K(0, 0);

  // extend list of integrals of monomials
  UpdateIntegrals_(c, 2 * order_);
  const Polynomial& integrals = integrals_[c];
   
  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt.end() <= p.end(); ++jt) {
      const int* idx_q = jt.multi_index();
      int l = jt.PolynomialPosition();
      
      int n(0);
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = idx_p[i] + idx_q[i];
        n += multi_index[i];
      }

      const auto& coefs = integrals.monomials(n).coefs();
      M(l, k) = M(k, l) = K00 * coefs[p.MonomialPosition(multi_index)];
    }
  }

  ChangeBasis_(c, M);

  return 0;
}


/* ******************************************************************
* Mass matrix for Taylor basis functions. 
****************************************************************** */
int DG_Modal::MassMatrix(
    int c, const Tensor& K, PolynomialOnMesh& integrals, DenseMatrix& M)
{
  double K00 = K(0, 0);

  // extend list of integrals of monomials
  numi_.UpdateMonomialIntegralsCell(c, 2 * order_, integrals);
   
  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt.end() <= p.end(); ++jt) {
      const int* idx_q = jt.multi_index();
      int l = jt.PolynomialPosition();
      
      int n(0);
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = idx_p[i] + idx_q[i];
        n += multi_index[i];
      }

      const auto& coefs = integrals.poly().monomials(n).coefs();
      M(l, k) = M(k, l) = K00 * coefs[p.MonomialPosition(multi_index)];
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
  Polynomial p(d_, order_);

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

      for (auto jt = it; jt.end() <= p.end(); ++jt) {
        const int* idx_q = jt.multi_index();
        int l = jt.PolynomialPosition();

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

  // symmetri part of mass matrix
  for (int k = 0; k < nrows; ++k) {
    for (int l = k + 1; l < nrows; ++l) {
      M(l, k) = M(k, l); 
    }
  }

  ChangeBasis_(c, M);

  return 0;
}


/* ******************************************************************
* Advection matrix for Taylor basis and cell-based velocity u.
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
  UpdateIntegrals_(c, order_ + std::max(0, order_ - 1) + uk);
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

    pgrad.Gradient(pp);
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

  // gradient operator is applied to solution
  if (!grad_on_test) {
    A.Transpose();
  }

  ChangeBasis_(c, A);

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
  mesh_->face_get_cells(f, (Parallel_type)WhetStone::USED, &cells);
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
      double vel1 = numi_.IntegratePolynomialsFace(f, polys);
      vel1 /= mesh_->face_area(f);
      vel1 *= dir;  

      // upwind-downwind integral
      polys[1] = &p1;

      double vel0 = numi_.IntegratePolynomialsFace(f, polys);
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

  // gradient operator is applied to solution
  if (!jump_on_test) {
    A.Transpose();
  }

  ChangeBasis_(cells[0], cells[1], A);

  return 0;
}


/* ******************************************************************
* Selects mean value of the polynomial defined on cell boundary.
* NOTE: This is a part of the elliptic projector, see MeshMaps.cc
****************************************************************** */
void DG_Modal::CoVelocityCell(
   int c, const std::vector<const Polynomial*> vf, VectorPolynomial& vc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  AmanziGeometry::Point u0(d_);
  WhetStone::Monomial mono(d_, 1);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    double area = mesh_->face_area(f);

    for (auto it = mono.begin(); it.end() <= 1; ++it) {
      Polynomial p(d_, it.multi_index());
      p.set_origin(xc);

      Polynomial q(*vf[n]);
      q *= 1.0 / area;
      q.ChangeOrigin(xc);
      p *= q;

      int k = it.MonomialPosition();
      u0[k] += dirs[n] * numi_.IntegratePolynomialFace(f, p);
    }
  }

  u0 /= mesh_->cell_volume(c);
  for (int i = 0; i < d_; ++i) {
    vc[i].ChangeOrigin(xc);
    vc[i](0, 0) = u0[i];
  }
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

  int nrows = A.NumRows();
  int m(nrows / 2);
  std::vector<double> a1(m), a2(m), b1(m), b2(m);

  Iterator it(d_);
  for (it.begin(); it.end() <= order_; ++it) {
    int k = it.PolynomialPosition();

    double ak, bk;
    TaylorBasis(c1, it, &ak, &bk);
    a1[k] = ak;
    b1[k] = -ak * bk;

    TaylorBasis(c2, it, &ak, &bk);
    a2[k] = ak;
    b2[k] = -ak * bk;
  }

  // calculate A * R
  for (int k = 1; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(i, k) = A(i, k) * a1[k] + A(i, 0) * b1[k];
      A(i, k + m) = A(i, k + m) * a2[k] + A(i, m) * b2[k];

      A(i + m, k) = A(i + m, k) * a1[k] + A(i + m, 0) * b1[k];
      A(i + m, k + m) = A(i + m, k + m) * a2[k] + A(i + m, m) * b2[k];
    }
  }

  // calculate R^T * A * R
  for (int k = 1; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(k, i) = A(k, i) * a1[k] + A(0, i) * b1[k];
      A(k + m, i) = A(k + m, i) * a2[k] + A(m, i) * b2[k];

      A(k, i + m) = A(k, i + m) * a1[k] + A(0, i + m) * b1[k];
      A(k + m, i + m) = A(k + m, i + m) * a2[k] + A(m, i + m) * b2[k];
    }
  }
}


/* ******************************************************************
* Transform monomial \psi_k to polynomial a (\psi_k - b \psi_0) where
* factor b orthogonolizes to constant and factor a normalizes to 1.
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
  if (order_ == 0) {
    Polynomial poly(d_, 0);
    poly(0, 0) = coefs[0];
    poly.set_origin(mesh_->cell_centroid(c));
    return poly;
  }

  ASSERT(scales_a_.size() != 0);
  ASSERT(scales_a_[c].size() == coefs.size());

  Polynomial poly(scales_a_[c]);
  poly.set_origin(mesh_->cell_centroid(c));

  auto jt = coefs.begin();
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    poly(m, k) *= *jt; 
    poly(0, 0) -= poly(m, k) * scales_b_[c](m, k);
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
    int ncells_wghost = mesh_->num_entities((Entity_kind)WhetStone::CELL,
                                            (Parallel_type)WhetStone::USED);
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
      numi_.IntegrateMonomialsCell(c, integrals_[c].monomials(k));
    }
  }
}


/* ******************************************************************
* Normalize and optionally orthogonalize Taylor basis functions.
****************************************************************** */
void DG_Modal::UpdateScales_(int c, int order)
{
  if (scales_a_.size() == 0) {
    int ncells_wghost = mesh_->num_entities((Entity_kind)WhetStone::CELL,
                                            (Parallel_type)WhetStone::USED);
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
          if (basis_ == TAYLOR_BASIS_NORMALIZED_ORTHO) {
            b = aux1[k] / volume;
          } else {
            b = 0.0;  // no orthogonalization to constants
          }

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

