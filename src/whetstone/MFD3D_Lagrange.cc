/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Lagrange-type element: degrees of freedom are ordered as follows:
    (1) nodal values in the natural order;
    (2) moments on faces groupped by face;
    (3) moments of edges, groupped by edge, if any;
    (4) moments inside cell.
*/

#include <cmath>
#include <vector>

// Amanzi
#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "Basis_Regularized.hh"
#include "CoordinateSystems.hh"
#include "GrammMatrix.hh"
#include "MFD3D_Lagrange.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Constructor parses the parameter list
 ****************************************************************** */
MFD3D_Lagrange::MFD3D_Lagrange(const Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : MFD3D(mesh), InnerProduct(mesh), use_lo_(false)
{
  order_ = plist.get<int>("method order");
  if (plist.isParameter("use low-order scheme"))
    use_lo_ = plist.get<bool>("use low-order scheme");
}


/* ******************************************************************
 * High-order consistency condition for the stiffness matrix.
 ****************************************************************** */
int
MFD3D_Lagrange::H1consistency2D_(int c, const Tensor& K, DenseMatrix& N,
                                 DenseMatrix& Ac)
{
  Kokkos::View<Entity_ID*> faces, nodes;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int nfaces = faces.extent(0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c, false);

  // calculate degrees of freedom
  Polynomial poly(d_, order_), pf, pc;
  if (order_ > 1) {
    pf.Reshape(d_ - 1, order_ - 2);
    pc.Reshape(d_, order_ - 2);
  }
  int nd = poly.size();
  int ndf = pf.size();
  int ndc = pc.size();

  int ndof = nnodes + nfaces * ndf + ndc;
  N.Reshape(ndof, nd);
  Ac.Reshape(ndof, ndof);

  R_.Reshape(ndof, nd);
  G_.Reshape(nd, nd);

  // pre-calculate integrals of monomials
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_ - 2, integrals_);

  // selecting regularized basis
  Basis_Regularized basis;
  basis.Init(mesh_, AmanziMesh::CELL, c, order_, integrals_.poly());

  // populate matrices N and R
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d_, index, factor);
    cmono.set_origin(xc);

    // N: degrees of freedom at vertices
    auto grad = Gradient(cmono);

    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(nnodes);

    AmanziGeometry::Point xv(d_);
    for (int i = 0; i < nnodes; i++) {
      int v = nodes(i);
      mesh_->node_get_coordinates(v, &xv);
      N(i, col) = cmono.Value(xv);
    }

    // N and R: degrees of freedom on faces
    for (int i = 0; i < nfaces; i++) {
      int f = faces(i);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      AmanziGeometry::Point normal = mesh_->face_normal(f);
      FaceCoordinateSystem(normal, tau);
      normal *= dirs(i);

      AmanziGeometry::Point conormal = K * normal;

      Kokkos::View<Entity_ID*> face_nodes;
      mesh_->face_get_nodes(f, face_nodes);
      int nfnodes = face_nodes.extent(0);

      if (order_ == 1 && col > 0) {
        for (int j = 0; j < nfnodes; j++) {
          int v = face_nodes(j);
          int pos = 0;
          for (pos = 0; pos < nodes.extent(0); ++pos) {
            if (nodes(pos) == v) { break; }
          }
          R_(pos, col) += factor * conormal[col - 1] / 2;
        }
      } else if (col > 0) {
        int v, pos0, pos1;
        AmanziGeometry::Point x0(d_), x1(d_), xm(d_), sm(d_);

        Polynomial tmp = grad * conormal;

        v = face_nodes(0);
        for (pos0 = 0; pos0 < nodes.extent(0); ++pos0) {
          if (nodes(pos0) == v) { break; }
        }
        mesh_->node_get_coordinates(v, &x0);

        v = face_nodes(1);
        for (pos1 = 0; pos1 < nodes.extent(0); ++pos1) {
          if (nodes(pos1) == v) { break; }
        }
        mesh_->node_get_coordinates(v, &x1);

        if (order_ == 2) {
          // Simpson rule with 3 points
          double q0 = tmp.Value(x0);
          double q1 = tmp.Value(x1);
          double qmid = tmp.Value(mesh_->face_centroid(f));

          R_(pos0, col) += (q0 - qmid) / 6;
          R_(pos1, col) += (q1 - qmid) / 6;
          R_(row, col) = qmid;
        } else if (order_ > 2) {
          if (col < 3) {
            // constant gradient contributes only to 0th moment
            R_(row, col) += tmp(0);
          } else {
            auto polys = ConvertMomentsToPolynomials_(order_);

            // Gauss-Legendre quadrature rule with (order_) points
            int m(order_ - 1);
            for (int n = 0; n < order_; ++n) {
              xm = x0 * q1d_points[m][n] + x1 * (1.0 - q1d_points[m][n]);
              sm[0] = 0.5 - q1d_points[m][n];

              double factor = q1d_weights[m][n] * tmp.Value(xm);
              R_(pos0, col) += polys[0].Value(sm) * factor;
              R_(pos1, col) += polys[1].Value(sm) * factor;

              for (int k = 0; k < m; ++k) {
                R_(row + k, col) += polys[k + 2].Value(sm) * factor;
              }
            }
          }
        }
      }

      if (order_ > 1) {
        for (auto jt = pf.begin(); jt < pf.end(); ++jt) {
          const int* jndex = jt.multi_index();
          Polynomial fmono(d_ - 1, jndex, 1.0);
          fmono.InverseChangeCoordinates(xf, tau);

          polys[1] = &fmono;

          int n = jt.PolynomialPosition();
          N(row + n, col) = numi.IntegratePolynomialsFace(f, polys) / area;
        }
        row += ndf;
      }
    }

    // N and R: degrees of freedom in cells
    if (cmono.order() > 1) {
      VectorPolynomial Kgrad = K * grad;
      Polynomial tmp = Divergence(Kgrad);

      for (auto jt = tmp.begin(); jt < tmp.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int n = jt.PolynomialPosition();

        R_(row + n, col) = -tmp(n) / basis.monomial_scales()[m] * volume;
      }
    }

    if (order_ > 1) {
      for (auto jt = pc.begin(); jt < pc.end(); ++jt) {
        int n = jt.PolynomialPosition();
        const int* jndex = jt.multi_index();

        int nm(0);
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        int m = MonomialSetPosition(d_, multi_index);
        double factor = basis.monomial_scales()[it.MonomialSetOrder()] *
                        basis.monomial_scales()[jt.MonomialSetOrder()];
        N(row + n, col) = integrals_.poly()(nm, m) * factor / volume;
      }
    }
  }

  // Gramm matrix for gradients of polynomials
  G_.Multiply(N, R_, true);

  // calculate R inv(G) R^T
  DenseMatrix RG(ndof, nd), Rtmp(nd, ndof);

  // to invert generate matrix, we add and subtruct positive number
  G_(0, 0) = 1.0;
  G_.Inverse();
  G_(0, 0) = 0.0;
  RG.Multiply(R_, G_, false);

  Rtmp.Transpose(R_);
  Ac.Multiply(RG, Rtmp, false);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * High-order consistency condition for the stiffness matrix.
 ****************************************************************** */
int
MFD3D_Lagrange::H1consistency3D_(int c, const Tensor& K, DenseMatrix& N,
                                 DenseMatrix& Ac)
{
  Kokkos::View<Entity_ID*> faces, edges, nodes, fnodes;
  std::vector<int> map;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  mesh_->cell_get_edges(c, edges);
  int nedges = edges.extent(0);

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int nfaces = faces.extent(0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c, false);

  // calculate degrees of freedom
  Polynomial poly(d_, order_), pf, pe, pc;
  if (order_ > 1) {
    pe.Reshape(d_ - 2, order_ - 2);
    pf.Reshape(d_ - 1, order_ - 2);
    pc.Reshape(d_, order_ - 2);
  }
  int nd = poly.size();
  int nde = pe.size();
  int ndf = pf.size();
  int ndc = pc.size();

  int ndof = nnodes + nedges * nde + nfaces * ndf + ndc;
  N.Reshape(ndof, nd);
  Ac.Reshape(ndof, ndof);

  R_.Reshape(ndof, nd);
  G_.Reshape(nd, nd);

  // pre-calculate integrals of monomials
  NumericalIntegration numi(mesh_);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_ - 2, integrals_);

  // selecting regularized basis
  Basis_Regularized basis;
  basis.Init(mesh_, AmanziMesh::CELL, c, order_, integrals_.poly());

  // define surface matrix
  // MFD3D_Lagrange mfd_surf(mesh_);

  // populate matrices N and R
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d_, index, factor);
    cmono.set_origin(xc);

    // N: degrees of freedom at vertices
    auto grad = Gradient(cmono);

    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(nnodes);

    AmanziGeometry::Point xv(d_);
    for (int i = 0; i < nnodes; i++) {
      int v = nodes(i);
      mesh_->node_get_coordinates(v, &xv);
      N(i, col) = cmono.Value(xv);
    }

    // N and R: degrees of freedom on faces and edges
    for (int i = 0; i < nfaces; i++) {
      int f = faces(i);
      AmanziGeometry::Point normal = mesh_->face_normal(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      FaceCoordinateSystem(normal, tau);
      normal *= dirs(i);
      AmanziGeometry::Point conormal = K * normal;

      Polynomial tmp = grad * conormal;
      tmp.ChangeCoordinates(xf, tau);

      // low-order moments are the degrees of freedom
      for (auto jt = pf.begin(); jt < pf.end(); ++jt) {
        const int* jndex = jt.multi_index();
        Polynomial fmono(d_ - 1, jndex, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);

        polys[1] = &fmono;

        int k = jt.PolynomialPosition();
        R_(row + k, col) = area * tmp(k);
        N(row + k, col) = numi.IntegratePolynomialsFace(f, polys) / area;
      }
      row += ndf;

      // the highest-order moments require H^1 projector on each face
      // -- we calculate some of its parts, bypassing polynomial machinery
      // mfd_surf.H1ConsistencySurface(f);
      // const DenseMatrix& R = mfd_surf.R();
      // const DenseMatrix& G = mfd_surf.G();

      // Add w^T * inv(G) * R^T to matrix R
      // Add moments to matrix N
    }

    // N and R: degrees of freedom in cells
    row += nfaces * ndf + nedges * nde;
    if (cmono.order() > 1) {
      VectorPolynomial Kgrad = K * grad;
      Polynomial tmp = Divergence(Kgrad);

      for (auto jt = tmp.begin(); jt < tmp.end(); ++jt) {
        int m = jt.MonomialSetOrder();
        int n = jt.PolynomialPosition();

        R_(row + n, col) = -tmp(n) / basis.monomial_scales()[m] * volume;
      }
    }

    if (order_ > 1) {
      for (auto jt = pc.begin(); jt < pc.end(); ++jt) {
        int n = jt.PolynomialPosition();
        const int* jndex = jt.multi_index();

        int nm(0);
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        int m = MonomialSetPosition(d_, multi_index);
        double factor = basis.monomial_scales()[it.MonomialSetOrder()] *
                        basis.monomial_scales()[jt.MonomialSetOrder()];
        N(row + n, col) = integrals_.poly()(nm, m) * factor / volume;
      }
    }
  }

  // Gramm matrix for gradients of polynomials
  G_.Multiply(N, R_, true);

  // calculate R inv(G) R^T
  DenseMatrix RG(ndof, nd), Rtmp(nd, ndof);

  // to invert generate matrix, we add and subtruct positive number
  G_(0, 0) = 1.0;
  G_.Inverse();
  G_(0, 0) = 0.0;
  RG.Multiply(R_, G_, false);

  Rtmp.Transpose(R_);
  Ac.Multiply(RG, Rtmp, false);

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * High-order consistency condition for the stiffness matrix.
 ****************************************************************** */
int
MFD3D_Lagrange::H1consistencySurface(int c, const Tensor& K, DenseMatrix& N,
                                     DenseMatrix& Ac)
{
  /*
    mesh_->face_to_cell_edge_map(f, c, &map);
    mesh_->face_get_edges_and_dirs(f, &fedges, &fdirs);
    int nfedges = fedges.size();

    for (int j = 0; j < nedges; ++j) {
      if (col < 3) {
        // constant gradient contributes only to 0th moment
        R_(row, col) += tmp(0);
      } else {
        int pos0, pos1;
        AmanziGeometry::Point x0(d_), x1(d_), xm(d_), sm(d_);

        auto polys = ConvertMomentsToPolynomials_(order_);

        // Gauss-Legendre quadrature rule with (order_) points
        int m(order_ - 1);
        for (int n = 0; n < order_; ++n) {
          xm = x0 * q1d_points[m][n] + x1 * (1.0 - q1d_points[m][n]);
          sm[0] = 0.5 - q1d_points[m][n];

          double factor = q1d_weights[m][n] * tmp.Value(xm);
          R_(pos0, col) += polys[0].Value(sm) * factor;
          R_(pos1, col) += polys[1].Value(sm) * factor;

          for (int k = 0; k < m; ++k) {
            R_(row + k, col) += polys[k + 2].Value(sm) * factor;
          }
        }
      }
    }
  */
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Stiffness matrix for a high-order scheme.
 ****************************************************************** */
int
MFD3D_Lagrange::StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Stiffness matrix on a manifold for a high-order scheme.
 ****************************************************************** */
int
MFD3D_Lagrange::StiffnessMatrixSurface(int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistencySurface(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
 * Generic projector on space of polynomials of order k in cell c.
 ****************************************************************** */
void
MFD3D_Lagrange::ProjectorCell_(int c, const std::vector<Polynomial>& vf,
                               const ProjectorType type,
                               const Polynomial* moments, Polynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  Kokkos::View<Entity_ID*> faces, nodes;

  mesh_->cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  mesh_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c, false);

  // calculate stiffness matrix.
  Tensor T(d_, 1);
  DenseMatrix A;

  T(0, 0) = 1.0;
  StiffnessMatrix(c, T, A);

  // number of degrees of freedom
  Polynomial pf;
  if (order_ > 1) pf.Reshape(d_ - 1, order_ - 2);

  int nd = PolynomialSpaceDimension(d_, order_);
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nnodes + nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  DenseVector vdof(ndof);
  std::vector<const PolynomialBase*> polys(2);
  NumericalIntegration numi(mesh_);

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mesh_, AmanziMesh::CELL, c, order_, ptmp);

  AmanziGeometry::Point xv(d_);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);

  int row(nnodes);

  // calculate DOFs on boundary
  for (int n = 0; n < nfaces; ++n) {
    int f = faces(n);

    Kokkos::View<Entity_ID*> face_nodes;
    mesh_->face_get_nodes(f, face_nodes);
    int nfnodes = face_nodes.extent(0);

    for (int j = 0; j < nfnodes; j++) {
      int v = face_nodes(j);
      mesh_->node_get_coordinates(v, &xv);

      int pos = 0;
      for (pos = 0; pos < nodes.extent(0); ++pos) {
        if (nodes(pos) == v) { break; }
      }
      vdof(pos) = vf[n].Value(xv);
    }

    if (order_ > 1) {
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      FaceCoordinateSystem(normal, tau);

      polys[0] = &(vf[n]);

      for (auto it = pf.begin(); it < pf.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }
  }

  // DOFs inside cell: copy moments from input data
  if (ndof_c > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    AMANZI_ASSERT(ndof_c == v3.NumRows());

    for (int n = 0; n < ndof_c; ++n) { vdof(row + n) = v3(n); }
  }

  // calculate polynomial coefficients (in vector v5)
  DenseVector v4(nd), v5(nd);
  R_.Multiply(vdof, v4, true);
  G_.Multiply(v4, v5, false);

  uc = basis.CalculatePolynomial(mesh_, c, order_, v5);

  // calculate the constant value for elliptic projector
  if (order_ == 1) {
    AmanziGeometry::Point grad(d_);
    for (int j = 0; j < d_; ++j) { grad[j] = uc(1, j); }

    double a1(0.0), a2(0.0), tmp;
    for (int n = 0; n < nfaces; ++n) {
      int f = faces(n);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);

      tmp = vf[n].Value(xf) - grad * (xf - xc);
      a1 += tmp * area;
      a2 += area;
    }

    uc(0) = a1 / a2;
  } else if (order_ >= 2) {
    v4 = integrals_.poly().coefs();
    basis.ChangeBasisMyToNatural(v4);
    v4.reshape(nd);
    uc(0) = vdof(row) - (v4 * v5) / volume;
  }

  // calculate L2 projector
  if (type == ProjectorType::L2 && ndof_c > 0) {
    v5(0) = uc(0);

    DenseMatrix M, M2;
    DenseVector v6(nd - ndof_c);
    Polynomial poly(d_, order_);
    NumericalIntegration numi(mesh_);

    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);
    GrammMatrix(poly, integrals_, basis, M);

    M2 = M.SubMatrix(ndof_c, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    const DenseVector& v3 = moments->coefs();
    for (int n = 0; n < ndof_c; ++n) {
      v4(n) = v3(n) * mesh_->cell_volume(c, false);
    }

    for (int n = 0; n < nd - ndof_c; ++n) { v4(ndof_c + n) = v6(n); }

    M.Inverse();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mesh_, c, order_, v5);
  }

  // set correct origin
  uc.set_origin(xc);
}


/* *****************************************************************
 * N and R are used. Here we use simplified versions.
 ***************************************************************** */
void
MFD3D_Lagrange::ProjectorCell_LO_(int c, const std::vector<Polynomial>& vf,
                                  Polynomial& uc)
{
  Kokkos::View<Entity_ID*> faces, nodes;
  Kokkos::View<int*> dirs;

  mesh_->cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  mesh_->cell_get_faces_and_dirs(c, faces, dirs);
  int num_faces = faces.extent(0);

  // populate matrix R (should be a separate routine lipnikov@lanl.gv)
  DenseMatrix R(nnodes, d_);
  AmanziGeometry::Point p(d_), pnext(d_), pprev(d_), v1(d_), v2(d_), v3(d_);

  R.PutScalar(0.0);
  for (int i = 0; i < num_faces; i++) {
    int f = faces(i);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    Kokkos::View<Entity_ID*> face_nodes;
    mesh_->face_get_nodes(f, face_nodes);
    int num_face_nodes = face_nodes.extent(0);

    for (int j = 0; j < num_face_nodes; j++) {
      int v = face_nodes[j];
      double u(0.5);

      if (d_ == 2) {
        u = 0.5 * dirs(i);
      } else {
        int jnext = (j + 1) % num_face_nodes;
        int jprev = (j + num_face_nodes - 1) % num_face_nodes;

        int vnext = face_nodes(jnext);
        int vprev = face_nodes(jprev);

        mesh_->node_get_coordinates(v, &p);
        mesh_->node_get_coordinates(vnext, &pnext);
        mesh_->node_get_coordinates(vprev, &pprev);

        v1 = pprev - pnext;
        v2 = p - fm;
        v3 = v1 ^ v2;
        u = dirs(i) * norm(v3) / (4 * area);
      }
      int pos = 0;
      for (pos = 0; pos < nodes.extent(0); ++pos) {
        if (nodes(pos) == v) { break; }
      }
      for (int k = 0; k < d_; k++) R(pos, k) += normal[k] * u;
    }
  }

  uc.Reshape(d_, 1, true);
  for (int i = 0; i < nnodes; i++) {
    for (int k = 0; k < d_; k++) { uc(k + 1) += R(i, k) * vf[i](0); }
  }
  uc *= 1.0 / mesh_->cell_volume(c, false);
}


/* ******************************************************************
 * Projector on the space of polynomials of order k in cell c.
 * Note: projector can be build only as a post-processor.
 ****************************************************************** */
void
MFD3D_Lagrange::ProjectorCellFromDOFs_(int c, const DenseVector& dofs,
                                       const ProjectorType type, Polynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  int nd = R_.NumCols();
  int ndof = R_.NumRows();
  AMANZI_ASSERT(ndof == dofs.NumRows() && nd > 0);

  double volume = mesh_->cell_volume(c, false);
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

  Kokkos::View<Entity_ID*> faces;
  mesh_->cell_get_faces(c, faces);
  int nfaces = faces.extent(0);
  int nnodes = nfaces;
  int ndof_c(ndof - nnodes);

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mesh_, AmanziMesh::CELL, c, order_, ptmp);

  // calculate polynomial coefficients (in vector v5)
  DenseVector v4(nd), v5(nd);
  R_.Multiply(dofs, v4, true);
  G_.Multiply(v4, v5, false);

  uc = basis.CalculatePolynomial(mesh_, c, order_, v5);

  // calculate the constant value for H1 projector
  if (order_ == 1) {
    AmanziGeometry::Point grad(d_);
    for (int j = 0; j < d_; ++j) grad[j] = uc(j + 1);

    double a1(0.0), a2(0.0), tmp;
    for (int n = 0; n < nfaces; ++n) {
      int f = faces(n);
      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

      int m = (n + 1) % nfaces;
      double tmp = (dofs(n) + dofs(m)) / 2 - grad * (xf - xc);
      a1 += tmp * area;
      a2 += area;
    }

    uc(0) = a1 / a2;
  } else if (order_ >= 2) {
    v4 = integrals_.poly().coefs();
    basis.ChangeBasisMyToNatural(v4);
    v4.reshape(nd);
    uc(0) = dofs(nnodes) - (v4 * v5) / volume;
  }

  // calculate L2 projector
  if (type == ProjectorType::L2 && ndof_c > 0) {
    v5(0) = uc(0);

    DenseMatrix M, M2;
    DenseVector v6(nd - ndof_c);
    Polynomial poly(d_, order_);
    NumericalIntegration numi(mesh_);

    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);
    GrammMatrix(poly, integrals_, basis, M);

    M2 = M.SubMatrix(ndof_c, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    for (int n = 0; n < ndof_c; ++n) { v4(n) = dofs(nnodes + n) * volume; }

    for (int n = 0; n < nd - ndof_c; ++n) { v4(ndof_c + n) = v6(n); }

    M.Inverse();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mesh_, c, order_, v5);
  }

  // set correct origin
  uc.set_origin(xc);
}


/* *****************************************************************
 * Convert basis (DOFs at end-points and moments) to basis of regular
 * polynomials on interval (-1/2, 1/2).
 ***************************************************************** */
std::vector<Polynomial>
MFD3D_Lagrange::ConvertMomentsToPolynomials_(int order)
{
  int n = order + 1;
  WhetStone::DenseMatrix T(n, n);
  T.PutScalar(0.0);

  // values at end points
  double b0, a0(1.0), a1(1.0);
  for (int i = 0; i < n; ++i) {
    T(0, i) = a0;
    T(1, i) = a1;
    a0 /= -2;
    a1 /= 2;
  }

  // moments of even power
  b0 = 1.0;
  for (int k = 2; k < n; k += 2) {
    a0 = b0;
    for (int i = 0; i < n; i += 2) {
      T(k, i) = a0 / (k + i - 1);
      a0 /= 4;
    }
    b0 /= 4;
  }

  // moments of odd power
  b0 = 0.25;
  for (int k = 3; k < n; k += 2) {
    a0 = b0;
    for (int i = 1; i < n; i += 2) {
      T(k, i) = a0 / (k + i - 1);
      a0 /= 4;
    }
    b0 /= 4;
  }

  T.Inverse();

  // convert columns of T to polynomials
  std::vector<Polynomial> polys(n);
  for (int k = 0; k < n; ++k) {
    polys[k].Reshape(1, order);
    for (int i = 0; i < n; ++i) polys[k](i) = T(i, k);
  }

  return polys;
}

} // namespace WhetStone
} // namespace Amanzi
