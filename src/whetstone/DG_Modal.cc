/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin methods.
*/

#include "Point.hh"

#include "BasisFactory.hh"
#include "DenseMatrix.hh"
#include "DG_Modal.hh"
#include "FunctionUpwind.hh"
#include "Monomial.hh"
#include "Polynomial.hh"
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneMeshUtils.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Constructor.
 ****************************************************************** */
DG_Modal::DG_Modal(const Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : numi_(mesh), mesh_(mesh), d_(mesh->space_dimension())
{
  order_ = plist.get<int>("method order");
  std::string basis_name = plist.get<std::string>("dg basis");

  int ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  basis_.resize(ncells_wghost);
  monomial_integrals_.resize(ncells_wghost);

  for (int c = 0; c < ncells_wghost; ++c) {
    monomial_integrals_[c].Reshape(d_, 0);
    monomial_integrals_[c](0) = mesh_->cell_volume(c, false);
  }

  BasisFactory factory;
  for (int c = 0; c < ncells_wghost; ++c) {
    basis_[c] = factory.Create(basis_name);
    basis_[c]->Init(mesh_, AmanziMesh::CELL, c, order_, monomial_integrals_[c]);
  }
}


/* ******************************************************************
 * Mass matrix for Taylor basis functions.
 ****************************************************************** */
int
DG_Modal::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
{
  double K00 = K(0, 0);

  // extend (optionally) the list of integrals of non-normalized monomials
  numi_.UpdateMonomialIntegralsCell(c, 2 * order_, monomial_integrals_[c]);

  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);

  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt < p.end(); ++jt) {
      const int* idx_q = jt.multi_index();
      int l = jt.PolynomialPosition();

      for (int i = 0; i < d_; ++i) { multi_index[i] = idx_p[i] + idx_q[i]; }

      M(l, k) = M(k, l) =
        K00 * monomial_integrals_[c](PolynomialPosition(d_, multi_index));
    }
  }

  basis_[c]->BilinearFormNaturalToMy(M);

  return 0;
}


/* ******************************************************************
 * Mass matrix for Taylor basis functions.
 ****************************************************************** */
int
DG_Modal::MassMatrix(int c, const Tensor& K, PolynomialOnMesh& integrals,
                     DenseMatrix& M)
{
  // tensor must be scalar for this call
  double K00 = K(0, 0);

  // extend list of integrals of monomials
  numi_.UpdateMonomialIntegralsCell(c, 2 * order_, integrals);

  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);

  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt < p.end(); ++jt) {
      const int* idx_q = jt.multi_index();
      int l = jt.PolynomialPosition();

      for (int i = 0; i < d_; ++i) { multi_index[i] = idx_p[i] + idx_q[i]; }

      M(k, l) = K00 * integrals.poly()(PolynomialPosition(d_, multi_index));
      M(l, k) = M(k, l);
    }
  }

  basis_[c]->BilinearFormNaturalToMy(M);

  return 0;
}


/* ******************************************************************
 * Mass matrix for Taylor basis and polynomial coefficient K.
 ****************************************************************** */
int
DG_Modal::MassMatrixPoly_(int c, const Polynomial& K, DenseMatrix& M)
{
  // rebase the polynomial
  Polynomial Kcopy(K);
  Kcopy.ChangeOrigin(mesh_->cell_centroid(c));

  // extend (optionally) the list of integrals of non-normalized monomials
  int uk(Kcopy.order());
  numi_.UpdateMonomialIntegralsCell(c, 2 * order_ + uk, monomial_integrals_[c]);

  // sum up integrals to the mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  M.Reshape(nrows, nrows);
  M.PutScalar(0.0);

  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto mt = Kcopy.begin(); mt < Kcopy.end(); ++mt) {
      const int* idx_K = mt.multi_index();
      int n = mt.PolynomialPosition();
      double factor = Kcopy(n);
      if (factor == 0.0) continue;

      for (auto jt = it; jt < p.end(); ++jt) {
        const int* idx_q = jt.multi_index();
        int l = jt.PolynomialPosition();

        for (int i = 0; i < d_; ++i) {
          multi_index[i] = idx_p[i] + idx_q[i] + idx_K[i];
        }

        M(k, l) +=
          factor * monomial_integrals_[c](PolynomialPosition(d_, multi_index));
      }
    }
  }

  // symmetric part of mass matrix
  for (int k = 0; k < nrows; ++k) {
    for (int l = k + 1; l < nrows; ++l) { M(l, k) = M(k, l); }
  }

  basis_[c]->BilinearFormNaturalToMy(M);

  return 0;
}


/* ******************************************************************
 * Mass matrix for Taylor basis and piecewise polynomial coefficient K.
 ****************************************************************** */
int
DG_Modal::MassMatrixPiecewisePoly_(int c, const VectorPolynomial& K,
                                   DenseMatrix& M)
{
  Kokkos::View<AmanziMesh::Entity_ID*> faces, nodes;
  mesh_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  // allocate memory for matrix
  Polynomial p(d_, order_);
  int nrows = p.size();
  M.Reshape(nrows, nrows);
  M.PutScalar(0.0);

  std::vector<const WhetStoneFunction*> polys(3);

  std::vector<AmanziGeometry::Point> xy(3);
  // xy[0] = cell_geometric_center(*mesh_, c);
  xy[0] = mesh_->cell_centroid(c);

  for (auto it = p.begin(); it < p.end(); ++it) {
    int k = it.PolynomialPosition();
    int s = it.MonomialSetOrder();
    const int* idx0 = it.multi_index();

    Monomial p0(d_, idx0, 1.0);
    p0.set_origin(xc);

    polys[0] = &p0;

    for (auto jt = it; jt < p.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = jt.PolynomialPosition();
      int t = jt.MonomialSetOrder();

      Monomial p1(d_, idx1, 1.0);
      p1.set_origin(xc);

      polys[1] = &p1;

      // sum up local contributions
      for (auto n = 0; n < nfaces; ++n) {
        int f = faces(n);
        mesh_->face_get_nodes(f, nodes);
        mesh_->node_get_coordinates(nodes(0), &(xy[1]));
        mesh_->node_get_coordinates(nodes(1), &(xy[2]));

        polys[2] = &(K[n]);

        M(k, l) +=
          numi_.IntegrateFunctionsSimplex(xy, polys, s + t + K[n].order());
      }
    }
  }

  // symmetric part of mass matrix
  for (int k = 0; k < nrows; ++k) {
    for (int l = k + 1; l < nrows; ++l) { M(l, k) = M(k, l); }
  }

  basis_[c]->BilinearFormNaturalToMy(M);

  return 0;
}


/* ******************************************************************
 * Stiffness matrix for Taylor basis functions.
 ****************************************************************** */
int
DG_Modal::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  Tensor Ktmp(d_, 2);
  if (K.rank() == 2) {
    Ktmp = K;
  } else {
    Ktmp.MakeDiagonal(K(0, 0));
  }

  // extend (optionally) the list of integrals of non-normalized monomials
  numi_.UpdateMonomialIntegralsCell(
    c, std::max(2 * order_ - 2, 0), monomial_integrals_[c]);

  // copy integrals to mass matrix
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  A.Reshape(nrows, nrows);

  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt < p.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();

      int multi_index[3];
      for (int i = 0; i < d_; ++i) { multi_index[i] = index[i] + jndex[i]; }

      double sum(0.0), tmp;
      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) {
          if (index[i] > 0 && jndex[j] > 0) {
            multi_index[i]--;
            multi_index[j]--;

            tmp = monomial_integrals_[c](PolynomialPosition(d_, multi_index));
            sum += Ktmp(i, j) * tmp * index[i] * jndex[j];

            multi_index[i]++;
            multi_index[j]++;
          }
        }
      }

      A(l, k) = A(k, l) = sum;
    }
  }

  basis_[c]->BilinearFormNaturalToMy(A);

  return 0;
}


/* ******************************************************************
 * Advection matrix for Taylor basis and cell-based velocity u.
 ****************************************************************** */
int
DG_Modal::AdvectionMatrixPoly_(int c, const VectorPolynomial& u, DenseMatrix& A,
                               bool grad_on_test)
{
  // rebase the polynomial
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  VectorPolynomial ucopy(u);
  ucopy.ChangeOrigin(xc);

  // extend (optionally) the list  of integrals of non-normalized monomials
  int order_tmp = order_ + std::max(order_ - 1, 0) + ucopy[0].order();
  numi_.UpdateMonomialIntegralsCell(c, order_tmp, monomial_integrals_[c]);

  // sum-up integrals to the advection matrix
  int multi_index[3];
  Polynomial p(d_, order_), q(d_, order_);

  int nrows = p.size();
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* idx_p = it.multi_index();
    int k = it.PolynomialPosition();

    // product of polynomials need to align origins
    Polynomial pp(d_, idx_p, 1.0);
    pp.set_origin(xc);

    Polynomial tmp(Gradient(pp) * ucopy);

    for (auto mt = tmp.begin(); mt < tmp.end(); ++mt) {
      const int* idx_K = mt.multi_index();
      int n = mt.PolynomialPosition();
      double factor = tmp(n);
      if (factor == 0.0) continue;

      for (auto jt = q.begin(); jt < q.end(); ++jt) {
        const int* idx_q = jt.multi_index();
        int l = PolynomialPosition(d_, idx_q);

        for (int i = 0; i < d_; ++i) { multi_index[i] = idx_q[i] + idx_K[i]; }

        A(k, l) +=
          factor * monomial_integrals_[c](PolynomialPosition(d_, multi_index));
      }
    }
  }

  // gradient operator is applied to solution
  if (!grad_on_test) { A.Transpose(); }

  basis_[c]->BilinearFormNaturalToMy(A);

  return 0;
}


/* ******************************************************************
 * Advection matrix for Taylor basis and piecewise polynomial velocity.
 ****************************************************************** */
int
DG_Modal::AdvectionMatrixPiecewisePoly_(int c, const VectorPolynomial& u,
                                        DenseMatrix& A, bool grad_on_test)
{
  Kokkos::View<Entity_ID*> faces, nodes;
  mesh_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  // rebase the velocity polynomial (due to dot-product)
  VectorPolynomial ucopy(u);
  ucopy.ChangeOrigin(xc);

  // allocate memory for matrix
  Polynomial p(d_, order_), q(d_, order_);

  int nrows = p.size();
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  std::vector<const WhetStoneFunction*> polys(2);

  std::vector<AmanziGeometry::Point> xy(3);
  // xy[0] = cell_geometric_center(*mesh_, c);
  xy[0] = mesh_->cell_centroid(c);

  for (auto it = p.begin(); it < p.end(); ++it) {
    int k = it.PolynomialPosition();
    int s = it.MonomialSetOrder();
    const int* idx0 = it.multi_index();

    Polynomial p0(d_, idx0, 1.0);
    p0.set_origin(xc);

    auto pgrad = Gradient(p0);

    for (auto jt = q.begin(); jt < q.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = jt.PolynomialPosition();
      int t = jt.MonomialSetOrder();

      Monomial p1(d_, idx1, 1.0);
      p1.set_origin(xc);

      polys[0] = &p1;

      // sum-up integrals to the advection matrix
      for (int n = 0; n < nfaces; ++n) {
        int f = faces(n);
        mesh_->face_get_nodes(f, nodes);
        mesh_->node_get_coordinates(nodes(0), &(xy[1]));
        mesh_->node_get_coordinates(nodes(1), &(xy[2]));

        Polynomial tmp(d_, 0);
        tmp.set_origin(xc);
        for (int i = 0; i < d_; ++i) { tmp += pgrad[i] * ucopy[n * d_ + i]; }
        polys[1] = &tmp;

        A(k, l) += numi_.IntegrateFunctionsSimplex(xy, polys, t + tmp.order());
      }
    }
  }

  // gradient operator is applied to solution
  if (!grad_on_test) { A.Transpose(); }

  basis_[c]->BilinearFormNaturalToMy(A);

  return 0;
}


/* ******************************************************************
 * Upwind/downwind matrix for Taylor basis and normal velocity u.n.
 * The normal in u.n is scaled by the face area. If jump_on_test=true,
 * we calculate
 *
 *   a \Int { (u.n) \rho^* [\psi] } dS
 *
 * where star means upwind/downwind, \psi is a test function, \rho is
 * a solution, a=-1 for upwind, a=1 for downwind. Otherwise, we
 * calculate
 *
 *   a \Int { (u.n) \psi^* [\rho] } dS
 ****************************************************************** */
int
DG_Modal::FluxMatrix(int f, const Polynomial& un, DenseMatrix& A, bool upwind,
                     bool jump_on_test)
{
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  int ncells = cells.extent(0);

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  // identify index of upwind/downwind cell (id)
  int dir, id(0);
  mesh_->face_normal(f, false, cells(0), &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  if (ncells > 1) {
    double vel = un.Value(xf) * dir;
    if (upwind) vel = -vel;

    if (vel > 0.0) {
      id = 1;
    } else {
      dir = -dir;
    }
  }

  int col(id * size);
  int row(size - col);

  // Calculate upwind/downwind cells
  int c1, c2;
  c1 = cells(id);

  if (ncells == 1) {
    c2 = c1;
  } else {
    c2 = cells(1 - id);
  }

  // integrate traces of polynomials on face f
  std::vector<const PolynomialBase*> polys(3);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Monomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    Monomial p1(d_, idx0, 1.0);
    p1.set_origin(mesh_->cell_centroid(c2));

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Monomial q(d_, idx1, 1.0);
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
  if (!jump_on_test) { A.Transpose(); }

  if (ncells == 1) {
    basis_[cells(0)]->BilinearFormNaturalToMy(A);
  } else {
    basis_[cells(0)]->BilinearFormNaturalToMy(
      basis_[cells(0)], basis_[cells(1)], A);
  }

  return 0;
}


/* ******************************************************************
 * Upwind/downwind at Gauss points to calculate matrix for the Taylor
 * basis and normal velocity u.n. The normal in u.n is scaled by the
 * face area. If jump_on_test=true, we calculate
 *
 *   a \Int { (u.n) \rho^* [\psi] } dS
 *
 * where star means upwind/downwind, \psi is a test function, \rho is
 * a solution, a=-1 for upwind, and a=1 for downwind. Otherwise, we
 * calculate
 *
 *   a \Int { (u.n) \psi^* [\rho] } dS
 ****************************************************************** */
int
DG_Modal::FluxMatrixGaussPoints(int f, const Polynomial& un, DenseMatrix& A,
                                bool upwind, bool jump_on_test)
{
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  int ncells = cells.extent(0);

  Polynomial poly(d_, order_);
  int size = poly.size();
  int order_tmp = 2 * order_ + un.order();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  // identify index of upwind/downwind cell (id)
  int dir;
  double area = mesh_->face_area(f);
  mesh_->face_normal(f, false, cells(0), &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  // Calculate integrals needed for scaling
  int c1, c2, pos0, pos1;
  if (!upwind) {
    dir = -dir;
    area = -area;
  }

  if (dir > 0) {
    c1 = cells(0);
    c2 = cells(ncells - 1);
    pos0 = 0;
  } else {
    c1 = cells(ncells - 1);
    c2 = cells(0);
    pos0 = size;
  }
  pos1 = size - pos0;

  // integrate traces of polynomials on face f
  std::vector<const WhetStoneFunction*> polys(3);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Monomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    Monomial p1(d_, idx0, 1.0);
    p1.set_origin(mesh_->cell_centroid(c2));

    FunctionUpwindPlus unplus(&un);
    FunctionUpwindMinus unminus(&un);

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Monomial q0(d_, idx1, 1.0);
      q0.set_origin(mesh_->cell_centroid(c1));

      Monomial q1(d_, idx1, 1.0);
      q1.set_origin(mesh_->cell_centroid(c2));

      double v00, v00p, v10m, v01p, v11m;
      if (ncells == 1) {
        polys[0] = &un;
        polys[1] = &p0;
        polys[2] = &q0;

        v00 =
          numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;
        A(k, l) = v00;
      } else {
        polys[0] = &unplus;
        polys[1] = &p0;
        polys[2] = &q0;
        v00p =
          numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        polys[0] = &unminus;
        polys[1] = &p1;
        v10m =
          numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        polys[2] = &q1;
        v11m =
          numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        polys[0] = &unplus;
        polys[1] = &p0;
        v01p =
          numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        A(pos0 + l, pos0 + k) = v00p;
        A(pos0 + l, pos1 + k) = v10m;
        A(pos1 + l, pos0 + k) = -v01p;
        A(pos1 + l, pos1 + k) = -v11m;
      }
    }
  }

  // jump operator is applied to solution
  if (!jump_on_test) { A.Transpose(); }

  if (ncells == 1) {
    basis_[cells(0)]->BilinearFormNaturalToMy(A);
  } else {
    basis_[cells(0)]->BilinearFormNaturalToMy(
      basis_[cells(0)], basis_[cells(1)], A);
  }

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
int
DG_Modal::FluxMatrixRusanov(int f, const VectorPolynomial& uc1,
                            const VectorPolynomial& uc2, const Polynomial& uf,
                            DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List nodes;
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, cells);
  int ncells = cells.extent(0);

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);
  if (ncells == 1) return 0; // FIXME

  // identify index of downwind cell (id)
  int dir;
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, cells(0), &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  // Calculate integrals needed for scaling
  int c1 = cells(0);
  int c2 = cells(1);

  // integrate traces of polynomials on face f
  normal *= -1;
  Polynomial uf1 = uc1 * normal;
  Polynomial uf2 = uc2 * normal;

  uf2.ChangeOrigin(uf1.origin());
  Polynomial ufn = (uf1 + uf2) * 0.5;

  double tmp = numi_.PolynomialMaxValue(f, ufn);
  tmp *= 0.5;
  uf1(0, 0) -= tmp;
  uf2(0, 0) += tmp;

  std::vector<const PolynomialBase*> polys(3);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Polynomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    Polynomial p1(d_, idx0, 1.0);
    p1.set_origin(mesh_->cell_centroid(c2));

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Polynomial q0(d_, idx1, 1.0);
      q0.set_origin(mesh_->cell_centroid(c1));

      Polynomial q1(d_, idx1, 1.0);
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

  basis_[cells(0)]->BilinearFormNaturalToMy(
    basis_[cells(0)], basis_[cells(1)], A);

  return 0;
}


/* *****************************************************************
 * Jump matrix for Taylor basis:
 *
 *   \Int_f ( {K \grad \rho} [\psi] ) dS
 ****************************************************************** */
int
DG_Modal::FaceMatrixJump(int f, const Tensor& K1, const Tensor& K2,
                         DenseMatrix& A)
{
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, cells);
  int ncells = cells.extent(0);

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells(0);
  int c2 = (ncells > 1) ? cells(1) : -1;

  // Calculate co-normals
  int dir;
  AmanziGeometry::Point conormal1(d_), conormal2(d_);
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, c1, &dir);

  normal /= norm(normal);
  conormal1 = K1 * normal;
  if (c2 >= 0) { conormal2 = K2 * normal; }

  // integrate traces of polynomials on face f
  double coef00, coef01, coef10, coef11;
  Polynomial p0, p1, q0, q1;
  VectorPolynomial pgrad;
  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Polynomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    pgrad = Gradient(p0);
    p0 = pgrad * conormal1;

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Polynomial q0(d_, idx1, 1.0);
      q0.set_origin(mesh_->cell_centroid(c1));

      polys[0] = &p0;
      polys[1] = &q0;
      coef00 = numi_.IntegratePolynomialsFace(f, polys);

      A(k, l) = coef00 / ncells;

      if (c2 >= 0) {
        Polynomial p1(d_, idx0, 1.0);
        p1.set_origin(mesh_->cell_centroid(c2));

        pgrad = Gradient(p1);
        p1 = pgrad * conormal2;

        Polynomial q1(d_, idx1, 1.0);
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
    basis_[c1]->BilinearFormNaturalToMy(A);
  } else {
    basis_[c1]->BilinearFormNaturalToMy(basis_[c1], basis_[c2], A);
  }

  return 0;
}


/* *****************************************************************
 * Penalty matrix for Taylor basis and penalty coefficient Kf
 * corresponding to the following integral:
 *
 *   \Int_f { K_f [\psi] [\rho] } dS
 ****************************************************************** */
int
DG_Modal::FaceMatrixPenalty(int f, double Kf, DenseMatrix& A)
{
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, cells);
  int ncells = cells.extent(0);

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells(0);
  int c2 = (ncells > 1) ? cells(1) : -1;

  // integrate traces of polynomials on face f
  double coef00, coef01, coef11;
  Polynomial p0, p1, q0, q1;
  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Polynomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Polynomial q0(d_, idx1, 1.0);
      q0.set_origin(mesh_->cell_centroid(c1));

      polys[0] = &p0;
      polys[1] = &q0;
      coef00 = numi_.IntegratePolynomialsFace(f, polys);

      A(k, l) = Kf * coef00;

      if (c2 >= 0) {
        Polynomial p1(d_, idx0, 1.0);
        p1.set_origin(mesh_->cell_centroid(c2));

        Polynomial q1(d_, idx1, 1.0);
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
    basis_[c1]->BilinearFormNaturalToMy(A);
  } else {
    basis_[c1]->BilinearFormNaturalToMy(basis_[c1], basis_[c2], A);
  }

  return 0;
}

} // namespace WhetStone
} // namespace Amanzi
