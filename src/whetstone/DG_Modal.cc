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

#include <tuple>

#include "Point.hh"

#include "BasisFactory.hh"
#include "DenseMatrix.hh"
#include "DG_Modal.hh"
#include "FunctionUpwind.hh"
#include "Monomial.hh"
#include "Polynomial.hh"
#include "SurfaceMeshLight.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneMeshUtils.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor.
****************************************************************** */
DG_Modal::DG_Modal(const Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
  : BilinearForm(mesh),
    numi_(mesh)
{
  order_ = plist.get<int>("method order");
  std::string basis_name = plist.get<std::string>("dg basis");

  numi_order_ = order_;
  if (plist.isParameter("quadrature order"))
    numi_order_ = plist.get<int>("quadrature order");

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  basis_.resize(ncells_wghost);
  monomial_integrals_.resize(ncells_wghost);

  for (int c = 0; c < ncells_wghost; ++c) {
    monomial_integrals_[c].Reshape(d_, 0);
    monomial_integrals_[c](0) = mesh_->cell_volume(c);
  }

  BasisFactory<AmanziMesh::MeshLight> factory;
  for (int c = 0; c < ncells_wghost; ++c) {
    basis_[c] = factory.Create(basis_name);
    basis_[c]->Init(mesh_, c, order_, monomial_integrals_[c]);
  }
}


/* ******************************************************************
* Schema.
****************************************************************** */
std::vector<SchemaItem> DG_Modal::schema() const
{
  int nk = PolynomialSpaceDimension(d_, order_);
  std::vector<SchemaItem> items;

  items.push_back(std::make_tuple(AmanziMesh::CELL, DOF_Type::MOMENT, nk));
  return items;
}


/* ******************************************************************
* Mass matrix for Taylor basis functions. 
****************************************************************** */
int DG_Modal::MassMatrix(int c, const Tensor& K, DenseMatrix& M)
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
      
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = idx_p[i] + idx_q[i];
      }

      M(l, k) = M(k, l) = K00 * monomial_integrals_[c](PolynomialPosition(d_, multi_index));
    }
  }

  basis_[c]->BilinearFormNaturalToMy(M);

  return 0;
}


/* ******************************************************************
* Mass matrix for Taylor basis functions. 
****************************************************************** */
int DG_Modal::MassMatrix(
    int c, const Tensor& K, PolynomialOnMesh& integrals, DenseMatrix& M)
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
      
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = idx_p[i] + idx_q[i];
      }

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
int DG_Modal::MassMatrix(int c, const Polynomial& K, DenseMatrix& M)
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

        M(k, l) += factor * monomial_integrals_[c](PolynomialPosition(d_, multi_index));
      }
    }
  }

  // symmetric part of mass matrix
  for (int k = 0; k < nrows; ++k) {
    for (int l = k + 1; l < nrows; ++l) {
      M(l, k) = M(k, l); 
    }
  }

  basis_[c]->BilinearFormNaturalToMy(M);

  return 0;
}


/* ******************************************************************
* Stiffness matrix for Taylor basis functions and tensorial coefficient
****************************************************************** */
int DG_Modal::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  Tensor Ktmp(d_, 2);
  if (K.rank() == 2) {
    Ktmp = K;
  } else {
    Ktmp.MakeDiagonal(K(0, 0));
  }

  // extend (optionally) the list of integrals of non-normalized monomials
  numi_.UpdateMonomialIntegralsCell(c, std::max(2 * order_ - 2, 0), monomial_integrals_[c]);

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
      
      for (int i = 0; i < d_; ++i) {
        multi_index[i] = index[i] + jndex[i];
      }

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
* Stiffness matrix for Taylor basis functions and polynomial coefficient
****************************************************************** */
int DG_Modal::StiffnessMatrix(int c, const MatrixPolynomial& K, DenseMatrix& A)
{
  int uk(0);
  for (int i = 0; i < d_; ++i) 
    for (int j = 0; j < d_; ++j) 
      uk = std::max(uk, K(i, j).order());

  // extend (optionally) the list of integrals of non-normalized monomials
  numi_.UpdateMonomialIntegralsCell(c, std::max(2 * order_ - 2 + uk, 0), monomial_integrals_[c]);

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
      
      double sum(0.0);
      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) {
          if (index[i] > 0 && jndex[j] > 0) {
            for (auto mt = K(i, j).begin(); mt < K(i, j).end(); ++mt) {
              const int* idx_K = mt.multi_index();
              int n = mt.PolynomialPosition();
              double factor = K(i, j)(n);

              for (int m = 0; m < d_; ++m) {
                multi_index[m] = index[m] + jndex[m] + idx_K[m];
              }

              multi_index[i]--;
              multi_index[j]--;

              factor *= monomial_integrals_[c](PolynomialPosition(d_, multi_index)); 
              sum += factor * index[i] * jndex[j];
            }
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
* Stiffness matrix for Taylor basis functions and general coefficient
****************************************************************** */
int DG_Modal::StiffnessMatrix(int c, const WhetStoneFunction* K, DenseMatrix& A)
{
  int multi_index[3];
  Polynomial p(d_, order_);

  int nrows = p.size();
  A.Reshape(nrows, nrows);

  std::vector<const WhetStoneFunction*> funcs(2);
  funcs[0] = K;

  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt < p.end(); ++jt) {
      const int* jndex = jt.multi_index();
      int l = jt.PolynomialPosition();

      for (int i = 0; i < d_; ++i) {
        multi_index[i] = index[i] + jndex[i];
      }

      double sum(0.0), tmp;
      for (int i = 0; i < d_; ++i) {
        if (index[i] > 0 && jndex[i] > 0) {
          multi_index[i] -= 2;

          Monomial p0(d_, multi_index, 1.0);
          p0.set_origin(mesh_->cell_centroid(c));
          funcs[1] = &p0;

          tmp = numi_.IntegrateFunctionsTriangulatedCell(c, funcs, numi_order_);
          sum += tmp * index[i] * jndex[i];

          multi_index[i] += 2;
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
int DG_Modal::AdvectionMatrix(
    int c, const VectorPolynomial& u, DenseMatrix& A, bool grad_on_test)
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

        for (int i = 0; i < d_; ++i) {
          multi_index[i] = idx_q[i] + idx_K[i];
        }

        A(k, l) += factor * monomial_integrals_[c](PolynomialPosition(d_, multi_index));
      }
    }
  }

  // gradient operator is applied to solution
  if (!grad_on_test) {
    A.Transpose();
  }

  basis_[c]->BilinearFormNaturalToMy(A);

  return 0;
}


/* ******************************************************************
* Upwind/downwind matrix for Taylor basis and normal velocity u.n. 
* The normal in u.n is scaled by the face area. If jump_on_test=true,
* we calculate
* 
*   \Int { (u.n) \rho^* [\psi] } dS
* 
* where star means upwind/downwind, \psi is a test function, \rho is 
* a solution. Otherwise, we calculate 
*
*   \Int { (u.n) \psi^* [\rho] } dS
****************************************************************** */
int DG_Modal::FluxMatrix(int f, const Polynomial& un, DenseMatrix& A,
                         bool upwind, bool jump_on_test, double* flux)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  int ncells = cells.size();

  int size = PolynomialSpaceDimension(d_, order_);
  PolynomialIterator it0(d_), it1(d_);
  it0.begin(0); 
  it1.begin(order_ + 1);

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  // identify index of upwind/downwind cell (id) 
  int dir, id(0); 
  mesh_->face_normal(f, false, cells[0], &dir);
  const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

  *flux = un.Value(xf) * dir;
  if (ncells > 1) {
    double vel = *flux;
    if (!upwind) vel = -vel;

    if (vel < 0.0) {
      id = 1;
    } else {
      dir = -dir;
    }
  }

  int col(id * size);
  int row(size - col);

  // Calculate upwind/downwind cells
  int c1, c2;
  c1 = cells[id];

  if (ncells == 1) {
    c2 = c1;
  } else {
    c2 = cells[1 - id];
  }

  // create integrator on a surface (used for 3D only)
  const AmanziGeometry::Point& normal = mesh_->face_normal(f);
  auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);

  Teuchos::RCP<const SurfaceMeshLight> surf_mesh = Teuchos::rcp(new SurfaceMeshLight(mesh_, f, *coordsys));
  NumericalIntegration numi_f(surf_mesh);

  // integrate traces of polynomials on face f
  std::vector<const PolynomialBase*> polys(3), polys0(2), polys1(2), polys_tmp(1);
  Polynomial un_tmp, p0_tmp, p1_tmp, q_tmp;

  if (d_ == 3) {
    polys_tmp[0] = &un;
    un_tmp = ConvertPolynomialsToSurfacePolynomial(xf, coordsys, polys_tmp);
  }

  for (auto it = it0; it < it1; ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    // add monomials to the product list
    Monomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    Monomial p1(d_, idx0, 1.0);
    p1.set_origin(mesh_->cell_centroid(c2));

    if (d_ == 3) {
      polys_tmp[0] = &p0;
      p0_tmp = ConvertPolynomialsToSurfacePolynomial(xf, coordsys, polys_tmp);
      p0_tmp *= un_tmp;
      polys0[0] = &p0_tmp;

      polys_tmp[0] = &p1;
      p1_tmp = ConvertPolynomialsToSurfacePolynomial(xf, coordsys, polys_tmp);
      p1_tmp *= un_tmp;
      polys1[0] = &p1_tmp;
    }

    for (auto jt = it0; jt < it1; ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Monomial q(d_, idx1, 1.0);
      q.set_origin(mesh_->cell_centroid(c1));

      // add monomial to the product list
      if (d_ == 3) {
        polys_tmp[0] = &q;
        q_tmp = ConvertPolynomialsToSurfacePolynomial(xf, coordsys, polys_tmp);
        polys0[1] = &q_tmp;
        polys1[1] = &q_tmp;
      }

      // downwind-downwind integral
      double vel0, vel1;
      if (d_ == 2) {
        polys[0] = &un;
        polys[1] = &p0;
        polys[2] = &q;
        vel1 = numi_.IntegratePolynomialsFace(f, polys);
      } else {
        vel1 = numi_f.IntegratePolynomialsCell(0, polys0);
      }
      vel1 /= mesh_->face_area(f);
      vel1 *= dir;  

      // upwind-downwind integral
      if (d_ == 2) {
        polys[1] = &p1;
        vel0 = numi_.IntegratePolynomialsFace(f, polys);
      } else {
        vel0 = numi_f.IntegratePolynomialsCell(0, polys1);
      }
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

  if (ncells == 1) {
    basis_[cells[0]]->BilinearFormNaturalToMy(A);
  } else { 
    basis_[cells[0]]->BilinearFormNaturalToMy(basis_[cells[0]], basis_[cells[1]], A);
  }

  return 0;
}


/* ******************************************************************
* Upwind/downwind at Gauss points to calculate matrix for the Taylor
* basis and normal velocity u.n. The normal in u.n is scaled by the 
* face area. If jump_on_test=true, we calculate
* 
*   \Int { (u.n) \rho^* [\psi] } dS
* 
* where star means upwind/downwind, \psi is a test function, \rho is 
* a solution. Otherwise, we calculate 
*
*   \Int { (u.n) \psi^* [\rho] } dS
****************************************************************** */
int DG_Modal::FluxMatrixGaussPoints(
    int f, const Polynomial& un, DenseMatrix& A, bool upwind, bool jump_on_test)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly(d_, order_);
  int size = poly.size();
  int order_tmp = 2 * order_ + un.order();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  // identify index of upwind/downwind cell (id) 
  int dir; 
  double area = mesh_->face_area(f);
  mesh_->face_normal(f, false, cells[0], &dir);

  // Calculate integrals needed for scaling
  int c1, c2, pos0, pos1;
  if (!upwind) {
    dir = -dir;
    area = -area;
  }

  if (dir > 0) {
    c1 = cells[0];
    c2 = cells[ncells - 1];
    pos0 = 0;
  } else {
    c1 = cells[ncells - 1];
    c2 = cells[0];
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

        v00 = numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;
        A(k, l) = v00;
      } else {
        polys[0] = &unplus;
        polys[1] = &p0;
        polys[2] = &q0;
        v00p = numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        polys[0] = &unminus;
        polys[1] = &p1;
        v10m = numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        polys[2] = &q1;
        v11m = numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        polys[0] = &unplus;
        polys[1] = &p0;
        v01p = numi_.IntegrateFunctionsTriangulatedFace(f, polys, order_tmp) / area;

        A(pos0 + l, pos0 + k) = v00p;
        A(pos0 + l, pos1 + k) = v10m;
        A(pos1 + l, pos0 + k) =-v01p;
        A(pos1 + l, pos1 + k) =-v11m;
      }
    }
  }

  // jump operator is applied to solution
  if (!jump_on_test) {
    A.Transpose();
  }

  if (ncells == 1) {
    basis_[cells[0]]->BilinearFormNaturalToMy(A);
  } else { 
    basis_[cells[0]]->BilinearFormNaturalToMy(basis_[cells[0]], basis_[cells[1]], A);
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
int DG_Modal::FluxMatrixRusanov(
    int f, const VectorPolynomial& uc1, const VectorPolynomial& uc2,
    const Polynomial& uf, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells, nodes;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);
  if (ncells == 1) return 0;  // FIXME

  // identify index of downwind cell (id)
  int dir; 
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, cells[0], &dir);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = cells[1];

  // integrate traces of polynomials on face f
  normal *= -1;
  Polynomial uf1 = uc1 * normal;
  Polynomial uf2 = uc2 * normal;

  uf2.ChangeOrigin(uf1.get_origin());
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

  basis_[cells[0]]->BilinearFormNaturalToMy(basis_[cells[0]], basis_[cells[1]], A);

  return 0;
}


/* *****************************************************************
* Jump matrix for Taylor basis using functions:
*   \Int_f ( {K \grad \rho} [\psi] ) dS
****************************************************************** */
int DG_Modal::FaceMatrixJump(
    int f, const WhetStoneFunction* K1, const WhetStoneFunction* K2, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = (ncells > 1) ? cells[1] : -1;

  // Calculate co-normals
  int dir;
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, c1, &dir);
  double area = mesh_->face_area(f);
  normal /= area;

  // integrate traces of polynomials on face f
  double coef00, coef01, coef10, coef11;
  VectorPolynomial pgrad;
  std::vector<const WhetStoneFunction*> funcs(3);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Polynomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    pgrad = Gradient(p0);
    p0 = pgrad * normal;

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Polynomial q0(d_, idx1, 1.0);
      q0.set_origin(mesh_->cell_centroid(c1));

      funcs[0] = &p0;
      funcs[1] = &q0;
      funcs[2] = K1;
      coef00 = numi_.IntegrateFunctionsTriangulatedFace(f, funcs, numi_order_);

      A(k, l) = coef00 / ncells;

      if (c2 >= 0) {
        Polynomial p1(d_, idx0, 1.0);
        p1.set_origin(mesh_->cell_centroid(c2));

        pgrad = Gradient(p1);
        p1 = pgrad * normal;

        Polynomial q1(d_, idx1, 1.0);
        q1.set_origin(mesh_->cell_centroid(c2));

        funcs[1] = &q1;
        coef01 = numi_.IntegrateFunctionsTriangulatedFace(f, funcs, numi_order_);

        funcs[0] = &p1;
        funcs[2] = K2;
        coef11 = numi_.IntegrateFunctionsTriangulatedFace(f, funcs, numi_order_);

        funcs[1] = &q0;
        coef10 = numi_.IntegrateFunctionsTriangulatedFace(f, funcs, numi_order_);

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
int DG_Modal::FaceMatrixPenalty(int f, double Kf, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = (ncells > 1) ? cells[1] : -1;

  // integrate traces of polynomials on face f
  double coef00, coef01, coef11;
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

}  // namespace WhetStone
}  // namespace Amanzi

