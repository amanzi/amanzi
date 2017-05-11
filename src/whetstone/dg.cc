/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin methods.
*/

#include "Point.hh"

#include "dg.hh"
#include "DenseMatrix.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Mass matrix for Taylor basis functions. 
****************************************************************** */
int DG::TaylorMassMatrix(int c, int order, DenseMatrix& M)
{
  // calculate monomials
  int n(1), m((2 * order + 1) * (order + 1));
  std::vector<int> shift(1, 0);
  DenseVector monomials(m);

  monomials.PutScalar(0.0);
  monomials(0) = mesh_->cell_volume(c);

  for (int k = 1; k <= 2 * order; ++k) {
    shift.resize(n, k - 1);
    IntegrateMonomialsCell_(c, k, &(monomials(n)));
    n += k + 1;
  }
   
  // copy integrals to mass matrix (2D algorithm)
  int nrows = M.NumRows();

  for (int k = 0; k < nrows; ++k) {
    for (int l = k; l < nrows; ++l) {
      int i = shift[k] * shift[l];
      M(k, l) = M(l, k) = monomials(k + l + i);
    }
  }

  return 0;
}


/* ******************************************************************
* Advection matrix for Taylor basis and polynomial velocity u.
****************************************************************** */
int DG::TaylorAdvectionMatrixCell(
    int c, int order, std::vector<AmanziGeometry::Point>& u, DenseMatrix& A)
{
  // calculate monomials
  int uk(std::pow(2 * u.size(), 0.5) - 1);
  int nk(1), mk((2 * (order + uk) + 1) * (order + uk));

  std::vector<int> pk(1, 0); 
  DenseVector monomials(mk + uk + 1);

  monomials.PutScalar(0.0);
  monomials(0) = mesh_->cell_volume(c);

  for (int k = 1; k < 2 * order + uk; ++k) {
    pk.push_back(nk);
    IntegrateMonomialsCell_(c, k, &(monomials(nk)));
    nk += k + 1;
  }
   
  // copy integrals to mass matrix
  int k1, l1, m1, mm, nrows(A.NumRows());

  A.PutScalar(0.0);

  // two loops for column polynomial
  for (int i = 0; i <= order; ++i) {
    for (int k = 0; k < i + 1; ++ k) {
      k1 = pk[i] + k;

      // two loops for coefficient
      for (int m = 0, m1 = 0; m <= uk; ++m) {
        for (int n = 0; n < m + 1; ++n, ++m1) {
          double ux(u[m1][0]), uy(u[m1][1]);

          // two loops for row polynomial
          for (int j = 1; j <= order; ++j) {
            mm = pk[i + m + j - 1] + k;
            l1 = pk[j];
            for (int l = 0; l < j + 1; ++l) {
              A(k1, l1 + l) += ux * monomials(mm + n + l) * (j - l);
            }

            for (int l = 1; l < j + 1; ++l) {
              A(k1, l1 + l) += uy * monomials(mm + n + l - 1) * l;
            }
          }
        }
      }
    }
  }

  return 0;
}


/* ******************************************************************
* Advection matrix for Taylor basis and constant velocity u (2D only)
* Velocity is given in the face-based Taylor basis.
****************************************************************** */
int DG::TaylorAdvectionMatrixFace(
    int f, int order, std::vector<AmanziGeometry::Point>& u, DenseMatrix& M)
{
  AmanziMesh::Entity_ID_List cells, nodes;

  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  M.PutScalar(0.0);

  // calculate normal velocity
  int dir, cd(cells[0]);
  std::vector<double> factors;
  const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cd, &dir);

  for (int i = 0; i < u.size(); ++i) {
    factors.push_back(u[i] * normal);
  }

  // identify downwind cell
  int id(0); 
  if (factors[0] > 0.0) {
    if (ncells == 1) return 0;
    cd = cells[1];
    id = 1;
  } else {
    for (auto it = factors.begin(); it != factors.end(); ++it) {
      *it *= -1.0;
    }
  }

  // integrate traces from downwind cell
  mesh_->face_get_nodes(f, &nodes);

  AmanziGeometry::Point x1(d_), x2(d_);
  mesh_->node_get_coordinates(nodes[0], &x1);
  mesh_->node_get_coordinates(nodes[1], &x2);

  const AmanziGeometry::Point& xc1 = mesh_->cell_centroid(cd);

  int nrows = M.NumRows() / 2;
  int ks(nrows * id), ls(nrows * (1 - id));
  int is(ks);

  for (int i = 0; i <= order; ++i) {
    for (int k = 0; k < i + 1; ++k) {
      int js(ks);
      for (int j = 0; j <= order; ++j) {
        for (int l = 0; l < j + 1; ++l) {
          M(is + k, js) += IntegrateMonomialsEdge_(x1, x2, i - k, k,
                                                           j - l, l, factors, xc1, xc1);
          js++;
        }
      }
    }
    is += i + 1;
  }

  // integrate traces from both cells
  if (ncells == 1) return 0;

  int c2(cells[0] + cells[1] - cd);
  const AmanziGeometry::Point& xc2 = mesh_->cell_centroid(c2);

  is = ks;
  for (int i = 0; i <= order; ++i) {
    for (int k = 0; k < i + 1; ++k) {
      int js(ls);
      for (int j = 0; j <= order; ++j) {
        for (int l = 0; l < j + 1; ++l) {
          M(is + k, js) -= IntegrateMonomialsEdge_(x1, x2, i - k, k,
                                                           j - l, l, factors, xc1, xc2);
          js++;
        }
      }
    }
    is += i + 1;
  }
}


/* ******************************************************************
* Integrate all monomials of order k in cell c.
* The cell must be star-shape w.r.t. to its centroid.
****************************************************************** */
void DG::IntegrateMonomialsCell_(int c, int k, double* monomials)
{
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
      IntegrateMonomialsFace_(f, k, tmp, monomials);
    } else if (d_ == 2) {
      mesh_->face_get_nodes(f, &nodes);

      AmanziGeometry::Point x1(d_), x2(d_);
      mesh_->node_get_coordinates(nodes[0], &x1);
      mesh_->node_get_coordinates(nodes[1], &x2);

      x1 -= xc;
      x2 -= xc;
      IntegrateMonomialsEdge_(x1, x2, k, tmp, monomials);
    }
  }
}


/* ******************************************************************
* Integrate all monomials of order k on face.
****************************************************************** */
void DG::IntegrateMonomialsFace_(int c, int k, double factor, double* monomials)
{
}


/* ******************************************************************
* Integrate all monomials of order k on edge via quadrature rules.
****************************************************************** */
void DG::IntegrateMonomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    int k, double factor, double* monomials)
{
  double a1; 
  AmanziGeometry::Point xm(d_);

  if (d_ == 2) {
    int m = k / 2;  // calculate quadrature rule

    for (int i = 0; i <= k; ++i) {
      for (int n = 0; n <= m; ++n) { 
        xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);
        a1 = std::pow(xm[0], k - i) * std::pow(xm[1], i);
        monomials[i] += factor * a1 * q1d_weights[m][n];      
      }
    }
  }
}


/* ******************************************************************
* Integrate two monomials of order k on edge via quadrature rules.
****************************************************************** */
double DG::IntegrateMonomialsEdge_(
    const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
    int ix, int iy, int jx, int jy, 
    const std::vector<double>& factors, 
    const AmanziGeometry::Point& xc1, const AmanziGeometry::Point& xc2)
{
  double a1, a2, a3, tmp(0.0); 
  AmanziGeometry::Point xm(d_), ym(d_), xe((x1 + x2) /2);

  if (d_ == 2) {
    int uk(std::pow(2 * factors.size(), 0.5) - 1);
    int m((ix + iy + jx + jy + uk) / 2);  // calculate quadrature rule

    for (int n = 0; n <= m; ++n) { 
      xm = x1 * q1d_points[m][n] + x2 * (1.0 - q1d_points[m][n]);

      a3 = factors[0];
      if (uk > 0) {  // FIXME
        ym = xm - xe;  
        a3 += factors[1] * ym[0] + factors[2] * ym[1];
      }      

      ym = xm - xc1;
      a1 = std::pow(ym[0], ix) * std::pow(ym[1], iy);

      xm -= xc2;
      a2 = std::pow(xm[0], jx) * std::pow(xm[1], jy);

      tmp += a1 * a2 * a3 * q1d_weights[m][n];      
    }
  }

  return tmp;
}


/* ******************************************************************
* Polynomial approximation of map x2 = F(x1).
* We assume that vectors of vertices have a proper length.
****************************************************************** */
int DG::TaylorLeastSquareFit(int order,
                             const std::vector<AmanziGeometry::Point>& x1, 
                             const std::vector<AmanziGeometry::Point>& x2,
                             std::vector<AmanziGeometry::Point>& u) const
{
  int nk = (order + 1) * (order + 2) / 2;
  int nx = x1.size();

  // evaluate basis functions at given points
  int i1(0);
  DenseMatrix psi(nk, nx);

  for (int k = 0; k <= order; ++k) {
    for (int i = 0; i < k + 1; ++i) {
      for (int n = 0; n < nx; ++n) {
        psi(i1, n) = std::pow(x1[n][0], k - i) * std::pow(x1[n][1], i);
      }
      i1++;
    }
  }
      
  // form linear system
  DenseMatrix A(nk, nk);
  DenseVector bx(nk), by(nk), ux(nk), uy(nk);

  for (int i = 0; i < nk; ++i) {
    for (int j = i; j < nk; ++j) {
      double tmp(0.0);
      for (int n = 0; n < nx; ++n) {
        tmp += psi(i, n) * psi(j, n);
      }
      A(i, j) = A(j, i) = tmp;
    }

    bx(i) = 0.0;
    by(i) = 0.0;
    for (int n = 0; n < nx; ++n) {
      bx(i) += x2[n][0] * psi(i, n);
      by(i) += x2[n][1] * psi(i, n);
    }
  }

  // solver linear systems
  A.Inverse();
  A.Multiply(bx, ux, false);
  A.Multiply(by, uy, false);

  u.clear();
  for (int i = 0; i < nk; ++i) {
    u.push_back(AmanziGeometry::Point(ux(i), uy(i)));
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

