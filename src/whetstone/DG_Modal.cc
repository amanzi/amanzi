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
  numi_.ChangeBasisRegularToNatural(c, Kcopy);

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

  // symmetric part of mass matrix
  for (int k = 0; k < nrows; ++k) {
    for (int l = k + 1; l < nrows; ++l) {
      M(l, k) = M(k, l); 
    }
  }

  ChangeBasis_(c, M);

  return 0;
}


/* ******************************************************************
* Stiffness matrix for Taylor basis functions. 
****************************************************************** */
int DG_Modal::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  Tensor Ktmp(d_, 2);
  if (K.rank() == 2) {
    Ktmp = K;
  } else {
    Ktmp.MakeDiagonal(K(0, 0));
  }

  double volume = mesh_->cell_volume(c);
  double scale = numi_.MonomialNaturalScale(1, volume);

  // extend list of integrals of monomials
  UpdateIntegrals_(c, 2 * order_ - 2);
  const Polynomial& integrals = integrals_[c];
   
  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  A.Reshape(nrows, nrows);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt.end() <= p.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();
      
      int n(0);
      int multi_index[3];
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = index[i] + jndex[i];
        n += multi_index[i];
      }

      double sum(0.0), tmp;
      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) {
          if (index[i] > 0 && jndex[j] > 0) {
            multi_index[i]--;
            multi_index[j]--;

            const auto& coefs = integrals.monomials(n - 2).coefs();
            tmp = coefs[p.MonomialPosition(multi_index)]; 
            sum += Ktmp(i, j) * tmp * index[i] * jndex[j];

            multi_index[i]++;
            multi_index[j]++;
          }
        }
      }

      A(l, k) = A(k, l) = sum * scale * scale; 
    }
  }

  ChangeBasis_(c, A);

  return 0;
}


/* ******************************************************************
* Advection matrix for Taylor basis and cell-based velocity u.
****************************************************************** */
int DG_Modal::AdvectionMatrixPoly(int c, const VectorPolynomial& u, DenseMatrix& A, bool grad_on_test)
{
  // rebase the polynomial
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  VectorPolynomial ucopy(u);
  for (int i = 0; i < d_; ++i) {
    ucopy[i].ChangeOrigin(xc);
    numi_.ChangeBasisRegularToNatural(c, ucopy[i]);
  }

  // extend list of integrals of monomials
  int uk(ucopy[0].order());
  UpdateIntegrals_(c, order_ + std::max(0, order_ - 1) + uk);
  const Polynomial& integrals = integrals_[c];

  // gradient of a naturally scaled polynomial needs correction
  double scale = numi_.MonomialNaturalScale(1, volume);

  // sum-up integrals to the advection matrix
  int multi_index[3];
  Polynomial p(d_, order_), q(d_, order_);
  VectorPolynomial pgrad;

  int nrows = p.size();
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    // product of polynomials need to align origins
    Polynomial pp(d_, idx_p, scale);
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
* Upwind flux matrix for Taylor basis and normal velocity u.n. 
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
int DG_Modal::FluxMatrixUpwind(int f, const Polynomial& un, DenseMatrix& A,
                               bool jump_on_test)
{
  AmanziMesh::Entity_ID_List cells, nodes;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
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
  int c1 = cells[id];
  int c2 = cells[1 - id];

  double volume1 = mesh_->cell_volume(c1);
  double volume2 = mesh_->cell_volume(c2);

  UpdateIntegrals_(c1, 2 * order_);
  UpdateIntegrals_(c2, 2 * order_);

  // integrate traces of polynomials on face f
  std::vector<const Polynomial*> polys(3);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);
    int s = it.MonomialOrder();

    double factor = numi_.MonomialNaturalScale(s, volume1);
    Polynomial p0(d_, idx0, factor);
    p0.set_origin(mesh_->cell_centroid(c1));

    factor = numi_.MonomialNaturalScale(s, volume2);
    Polynomial p1(d_, idx0, factor);
    p1.set_origin(mesh_->cell_centroid(c2));

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);
      int t = jt.MonomialOrder();

      factor = numi_.MonomialNaturalScale(t, volume1);
      Polynomial q(d_, idx1, factor);
      q.set_origin(mesh_->cell_centroid(c1));

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

  // jump operator is applied to solution
  if (!jump_on_test) {
    A.Transpose();
  }

  ChangeBasis_(cells[0], cells[1], A);

  return 0;
}


/* ******************************************************************
* Rusanov flux matrix for Taylor basis and normal velocity u.n. 
* Velocities are given in the face-based Taylor basis. We calculate
* 
*   \Int { (u.n \rho)^* [\psi] } dS
*
* where (u.n \rho)^* is the Rusanov flux.
****************************************************************** */
int DG_Modal::FluxMatrixRusanov(
    int f, const VectorPolynomial& uc1, const VectorPolynomial& uc2,
    const Polynomial& uf, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells, nodes;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly0(d_, order_), poly1(d_, order_);
  int size = poly0.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);
  if (ncells == 1) return 0;  // FIXME

  // identify index of downwind cell (id)
  int dir; 
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, cells[0], &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = cells[1];

  double volume1 = mesh_->cell_volume(c1);
  double volume2 = mesh_->cell_volume(c2);

  UpdateIntegrals_(c1, 2 * order_);
  UpdateIntegrals_(c2, 2 * order_);

  // integrate traces of polynomials on face f
  normal *= -1;
  Polynomial uf1 = uc1 * normal;
  Polynomial uf2 = uc2 * normal;

  Polynomial ufn = (uf1 + uf2) * 0.5;

  double tmp = numi_.PolynomialMaxValue(f, ufn);
  tmp *= 0.5;
  uf1(0, 0) -= tmp;
  uf2(0, 0) += tmp;

  std::vector<const Polynomial*> polys(3);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);
    int s = it.MonomialOrder();

    double factor = numi_.MonomialNaturalScale(s, volume1);
    Polynomial p0(d_, idx0, factor);
    p0.set_origin(mesh_->cell_centroid(c1));

    factor = numi_.MonomialNaturalScale(s, volume2);
    Polynomial p1(d_, idx0, factor);
    p1.set_origin(mesh_->cell_centroid(c2));

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);
      int t = jt.MonomialOrder();

      factor = numi_.MonomialNaturalScale(t, volume1);
      Polynomial q0(d_, idx1, factor);
      q0.set_origin(mesh_->cell_centroid(c1));

      factor = numi_.MonomialNaturalScale(t, volume2);
      Polynomial q1(d_, idx1, factor);
      q1.set_origin(mesh_->cell_centroid(c2));

      double coef00, coef01, coef10, coef11;
      double scale = 2 * mesh_->face_area(f);

      // upwind-upwind integral
      polys[0] = &uf1;
      polys[1] = &p0;
      polys[2] = &q0;
      coef00 = numi_.IntegratePolynomialsFace(f, polys);

      // upwind-downwind integral
      polys[2] = &q1;
      coef01 = numi_.IntegratePolynomialsFace(f, polys);

      // downwind-downwind integral
      polys[0] = &uf2;
      polys[1] = &p1;
      coef11 = numi_.IntegratePolynomialsFace(f, polys);

      // downwind-upwind integral
      polys[2] = &q0;
      coef10 = numi_.IntegratePolynomialsFace(f, polys);

      A(l, k) = coef00 / scale;
      A(size + l, k) = -coef01 / scale;
      A(l, size + k) = coef10 / scale;
      A(size + l, size + k) = -coef11 / scale;
    }
  }

  ChangeBasis_(cells[0], cells[1], A);

  return 0;
}


/* *****************************************************************
* Jump matrix for Taylor basis:
*
*   \Int_f ( {K \grad \rho} [\psi] ) dS
****************************************************************** */
int DG_Modal::FaceMatrixJump(int f, const Tensor& K1, const Tensor& K2, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly0(d_, order_), poly1(d_, order_);
  int size = poly0.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = (ncells > 1) ? cells[1] : -1;

  double volume1 = mesh_->cell_volume(c1);
  double volume2(0.0);

  UpdateIntegrals_(c1, 2 * order_ - 1);
  if (c2 >= 0) {
    volume2 = mesh_->cell_volume(c2);
    UpdateIntegrals_(c2, 2 * order_ - 1);
  }

  // Calculate co-normals
  int dir;
  AmanziGeometry::Point conormal1(d_), conormal2(d_);
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, c1, &dir);

  normal /= norm(normal);
  conormal1 = K1 * normal;
  if (c2 >= 0) {
    conormal2 = K2 * normal;
  }

  // integrate traces of polynomials on face f
  double coef00, coef01, coef10, coef11;
  Polynomial p0, p1, q0, q1;
  VectorPolynomial pgrad;
  std::vector<const Polynomial*> polys(2);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);
    int s = it.MonomialOrder();

    double factor = numi_.MonomialNaturalScale(s, volume1);
    Polynomial p0(d_, idx0, factor);
    p0.set_origin(mesh_->cell_centroid(c1));

    pgrad.Gradient(p0);
    p0 = pgrad * conormal1;

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);
      int t = jt.MonomialOrder();

      factor = numi_.MonomialNaturalScale(t, volume1);
      Polynomial q0(d_, idx1, factor);
      q0.set_origin(mesh_->cell_centroid(c1));

      polys[0] = &p0;
      polys[1] = &q0;
      coef00 = numi_.IntegratePolynomialsFace(f, polys);

      A(k, l) = coef00 / ncells;

      if (c2 >= 0) {
        factor = numi_.MonomialNaturalScale(s, volume2);
        Polynomial p1(d_, idx0, factor);
        p1.set_origin(mesh_->cell_centroid(c2));

        pgrad.Gradient(p1);
        p1 = pgrad * conormal2;

        factor = numi_.MonomialNaturalScale(t, volume2);
        Polynomial q1(d_, idx1, factor);
        q1.set_origin(mesh_->cell_centroid(c2));

        polys[1] = &q1;
        coef01 = numi_.IntegratePolynomialsFace(f, polys);

        polys[0] = &p1;
        coef11 = numi_.IntegratePolynomialsFace(f, polys);

        polys[1] = &q0;
        coef10 = numi_.IntegratePolynomialsFace(f, polys);

        A(k, size + l) = -coef01 / ncells;
        A(size + k, size + l) = -coef11 / ncells;
        A(size + k, l) = coef10 / ncells;
      }
    }
  }

  if (ncells == 1) {
    ChangeBasis_(c1, A);
  } else {
    ChangeBasis_(c1, c2, A);
  }
}


/* *****************************************************************
* Penalty matrix for Taylor basis and penalty coefficient Kf
* corresponding to the following integral:
*
*   \Int_f { K_f [\psi] [\rho] } dS
****************************************************************** */
int DG_Modal::FaceMatrixPenalty(int f, double Kf, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly0(d_, order_), poly1(d_, order_);
  int size = poly0.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = (ncells > 1) ? cells[1] : -1;

  double volume1 = mesh_->cell_volume(c1);
  double volume2(0.0);

  UpdateIntegrals_(c1, 2 * order_);
  if (c2 >= 0) {
    volume2 = mesh_->cell_volume(c2);
    UpdateIntegrals_(c2, 2 * order_);
  }

  // integrate traces of polynomials on face f
  double coef00, coef01, coef11;
  Polynomial p0, p1, q0, q1;
  std::vector<const Polynomial*> polys(2);

  for (auto it = poly0.begin(); it.end() <= poly0.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = poly0.PolynomialPosition(idx0);
    int s = it.MonomialOrder();

    double factor = numi_.MonomialNaturalScale(s, volume1);
    Polynomial p0(d_, idx0, factor);
    p0.set_origin(mesh_->cell_centroid(c1));

    for (auto jt = poly1.begin(); jt.end() <= poly1.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = poly1.PolynomialPosition(idx1);
      int t = jt.MonomialOrder();

      factor = numi_.MonomialNaturalScale(t, volume1);
      Polynomial q0(d_, idx1, factor);
      q0.set_origin(mesh_->cell_centroid(c1));

      polys[0] = &p0;
      polys[1] = &q0;
      coef00 = numi_.IntegratePolynomialsFace(f, polys);

      A(k, l) = Kf * coef00;

      if (c2 >= 0) {
        factor = numi_.MonomialNaturalScale(s, volume2);
        Polynomial p1(d_, idx0, factor);
        p1.set_origin(mesh_->cell_centroid(c2));

        factor = numi_.MonomialNaturalScale(t, volume2);
        Polynomial q1(d_, idx1, factor);
        q1.set_origin(mesh_->cell_centroid(c2));

        polys[1] = &q1;
        coef01 = numi_.IntegratePolynomialsFace(f, polys);

        polys[0] = &p1;
        coef11 = numi_.IntegratePolynomialsFace(f, polys);

        A(k, size + l) = -Kf * coef01;
        A(size + k, size + l) = Kf * coef11;
        A(size + l, k) = -Kf * coef01;
      }
    }
  }

  if (ncells == 1) {
    ChangeBasis_(c1, A);
  } else {
    ChangeBasis_(c1, c2, A);
  }

  return 0;
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
  } else if (basis_ == TAYLOR_BASIS_NATURAL) {
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
Polynomial DG_Modal::CalculatePolynomial(int c, const DenseVector& coefs) const
{
  if (order_ == 0 || basis_ == TAYLOR_BASIS_NATURAL) {
    Polynomial poly(d_, order_);
    poly.SetPolynomialCoefficients(coefs);
    poly.set_origin(mesh_->cell_centroid(c));
    return poly;
  }

  ASSERT(scales_a_.size() != 0);
  ASSERT(scales_a_[c].size() == coefs.NumRows());

  Polynomial poly(scales_a_[c]);
  poly.set_origin(mesh_->cell_centroid(c));

  int n(0);
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    poly(m, k) *= coefs(n++); 
    poly(0, 0) -= poly(m, k) * scales_b_[c](m, k);
  }

  return poly;
}


/* ******************************************************************
* Update integrals of non-normalized monomials.
****************************************************************** */
void DG_Modal::UpdateIntegrals_(int c, int order)
{
  if (integrals_.size() == 0) {
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
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
    int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
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

