/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Lagrange-type element: degrees of freedom are nodal values and
  moments on edges, faces and inside cell.
*/

#include <cmath>
#include <vector>

// Amanzi
#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

// WhetStone
#include "CoordinateSystems.hh"
#include "DG_Modal.hh"
#include "MFD3D_Lagrange.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
* Only the upper triangular part of Ac is calculated. 
****************************************************************** */
int MFD3D_Lagrange::H1consistency(
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  ASSERT(d_ == 2);  // FIXME
  ASSERT(order_ < 4);

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c); 
  double volume = mesh_->cell_volume(c); 

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

  // populate matrices N and R
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const Polynomial*> polys(2);

  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) { 
    const int* index = it.multi_index();
    double factor = numi.MonomialNaturalScale(it.MonomialOrder(), volume);
    Polynomial cmono(d_, index, factor);
    cmono.set_origin(xc);  

    // N and R: degrees of freedom on faces 
    VectorPolynomial grad;
    grad.Gradient(cmono);
     
    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(nnodes);

    AmanziGeometry::Point xv(d_);
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      mesh_->node_get_coordinates(v, &xv);
      N(i, col) = cmono.Value(xv);
    }

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      AmanziGeometry::Point normal = mesh_->face_normal(f);
      FaceCoordinateSystem(normal, tau);
      normal *= dirs[i];

      Entity_ID_List face_nodes;
      mesh_->face_get_nodes(f, &face_nodes);
      int nfnodes = face_nodes.size();

      if (order_ == 1 && col > 0) {
        for (int j = 0; j < nfnodes; j++) {
          int v = face_nodes[j];
          int pos = FindPosition(v, nodes);
          R_(pos, col) += factor * normal[col - 1] / 2;
        }
      } else if (col > 0) {
        int v, pos0, pos1;
        AmanziGeometry::Point x0(d_), x1(d_), xm(d_), sm(d_);

        Polynomial tmp = grad * normal;

        v = face_nodes[0];
        pos0 = FindPosition(v, nodes);
        mesh_->node_get_coordinates(v, &x0);

        v = face_nodes[1];
        pos1 = FindPosition(v, nodes);
        mesh_->node_get_coordinates(v, &x1);

        if (order_ == 2) {
          // Simpson rule with 3 points
          double q0 = tmp.Value(x0);
          double q1 = tmp.Value(x1);
          double qmid = tmp.Value(mesh_->face_centroid(f));

          R_(pos0, col) += (q0 - qmid) / 6;
          R_(pos1, col) += (q1 - qmid) / 6;
          R_(row,  col) = qmid;
        } else if (order_ == 3) {
          // Gauss-Legendre quadrature rule with 3 points
          int m(2); 
          Polynomial poly0(1, 3), poly1(1, 3), poly2(1, 2), poly3(1, 3);

          poly0(0, 0) = -0.25;
          poly0(1, 0) = 1.5;
          poly0(2, 0) = 3.0;
          poly0(3, 0) = -10.0;

          poly1(0, 0) = -0.25;
          poly1(1, 0) = -1.5;
          poly1(2, 0) = 3.0;
          poly1(3, 0) = 10.0;

          poly2(0, 0) = 1.5;
          poly2(1, 0) = 0.0;
          poly2(2, 0) = -6.0;

          poly3(0, 0) = 0.0;
          poly3(1, 0) = 30.0;
          poly3(2, 0) = 0.0;
          poly3(3, 0) =-120.0;

          for (int n = 0; n < 3; ++n) { 
            xm = x0 * q1d_points[m][n] + x1 * (1.0 - q1d_points[m][n]);
            sm[0] = 0.5 - q1d_points[m][n];

            double factor = q1d_weights[m][n] * tmp.Value(xm);
            R_(pos0, col) += poly0.Value(sm) * factor;
            R_(pos1, col) += poly1.Value(sm) * factor;

            R_(row, col) += poly2.Value(sm) * factor;
            R_(row + 1, col) += poly3.Value(sm) * factor;
          }
        }
      }

      if (order_ > 1) {
        for (auto jt = pf.begin(); jt.end() <= pf.end(); ++jt) {
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
      Polynomial tmp = cmono.Laplacian();
      numi.ChangeBasisRegularToNatural(c, tmp);

      for (auto jt = tmp.begin(); jt.end() <= tmp.end(); ++jt) {
        int m = jt.MonomialOrder();
        int k = jt.MonomialPosition();
        int n = jt.PolynomialPosition();

        R_(row + n, col) = -tmp(m, k) * volume;
      }
    }

    if (order_ > 1) {
      for (auto jt = pc.begin(); jt.end() <= pc.end(); ++jt) {
        int n = jt.PolynomialPosition();
        const int* jndex = jt.multi_index();

        int nm(0);
        int multi_index[3];
        for (int i = 0; i < d_; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        const auto& coefs = integrals_.poly().monomials(nm).coefs();
        N(row + n, col) = coefs[poly.MonomialPosition(multi_index)] / volume; 
      }
    }
  }

  // set the Gramm-Schidt matrix for gradients of polynomials
  G_.PutScalar(0.0);

  // -- gradient of a naturally scaled polynomial needs correction
  double scale = numi.MonomialNaturalScale(1, volume);
   
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    const int* index = it.multi_index();
    int k = it.PolynomialPosition();

    for (auto jt = it; jt.end() <= poly.end(); ++jt) {
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
        if (index[i] > 0 && jndex[i] > 0) {
          multi_index[i] -= 2;
          const auto& coefs = integrals_.poly().monomials(n - 2).coefs();
          tmp = coefs[poly.MonomialPosition(multi_index)]; 
          sum += tmp * index[i] * jndex[i];
          multi_index[i] += 2;
        }
      }

      G_(l, k) = G_(k, l) = K(0, 0) * sum * scale * scale; 
    }
  }

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
* Stiffness matrix for a high-order scheme.
****************************************************************** */
int MFD3D_Lagrange::StiffnessMatrix(
    int c, const Tensor& K, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, K, N, A);
  if (ok) return ok;

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
void MFD3D_Lagrange::ProjectorCell_(
    int c, const std::vector<VectorPolynomial>& vf, 
    const Projectors::Type type, bool is_harmonic,
    const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) 
{
  ASSERT(d_ == 2);

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // calculate stiffness matrix.
  Tensor T(d_, 1);
  DenseMatrix A;

  T(0, 0) = 1.0;
  StiffnessMatrix(c, T, A);  

  // number of degrees of freedom
  Polynomial pf;
  if (order_ > 1) {
    pf.Reshape(d_ - 1, order_ - 2);
  }

  int nd = G_.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nnodes + nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  // auxiliary data for optional calculation of internal moments
  DenseMatrix Acf, Acc;
  if (ndof_c > 0 && is_harmonic) {
    Acf = A.SubMatrix(ndof_f, ndof, 0, ndof_f);
    Acc = A.SubMatrix(ndof_f, ndof, ndof_f, ndof);
    Acc.Inverse();
  }
  
  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, order_, true);
  }

  DenseVector vdof(ndof);
  std::vector<const Polynomial*> polys(2);
  NumericalIntegration numi(mesh_);

  AmanziGeometry::Point xv(d_);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);

  for (int i = 0; i < dim; ++i) {
    int row(nnodes);

    // calculate DOFs on boundary
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];

      Entity_ID_List face_nodes;
      mesh_->face_get_nodes(f, &face_nodes);
      int nfnodes = face_nodes.size();

      for (int j = 0; j < nfnodes; j++) {
        int v = face_nodes[j];
        mesh_->node_get_coordinates(v, &xv);

        int pos = WhetStone::FindPosition(v, nodes);
        vdof(pos) = vf[n][i].Value(xv);
      }

      if (order_ > 1) { 
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
        double area = mesh_->face_area(f);

        // local coordinate system with origin at face centroid
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        FaceCoordinateSystem(normal, tau);

        polys[0] = &(vf[n][i]);

        for (auto it = pf.begin(); it.end() <= pf.end(); ++it) {
          const int* index = it.multi_index();
          Polynomial fmono(d_ - 1, index, 1.0);
          fmono.InverseChangeCoordinates(xf, tau);  

          polys[1] = &fmono;

          vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
          row++;
        }
      }
    }

    // calculate DOFs inside cell using
    // -- either harmonic extension inside cell
    if (ndof_c > 0 && is_harmonic) {
      DenseVector v1(ndof_f), v2(ndof_c), v3(ndof_c);

      for (int n = 0; n < ndof_f; ++n) {
        v1(n) = vdof(n);
      }

      Acf.Multiply(v1, v2, false);
      Acc.Multiply(v2, v3, false);

      moments->Reshape(ndof_c);
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = -v3(n);
        (*moments)(n) = -v3(n);
      }
    }
    // -- or copy moments from input data
    else if (ndof_c > 0 && !is_harmonic) {
      ASSERT(ndof_c == moments->NumRows());
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = (*moments)(n);
      }
    }

    // calculate polynomial coefficients (in vector v5)
    DenseVector v4(nd), v5(nd);
    R_.Multiply(vdof, v4, true);
    G_.Multiply(v4, v5, false);

    uc[i].SetPolynomialCoefficients(v5);
    numi.ChangeBasisNaturalToRegular(c, uc[i]);

    // calculate the constant value for elliptic projector
    if (order_ == 1) {
      AmanziGeometry::Point grad(d_), zero(d_);
      for (int j = 0; j < d_; ++j) {
        grad[j] = uc[i](1, j);
      }
    
      double a1(0.0), a2(0.0), tmp;
      for (int n = 0; n < nfaces; ++n) {  
        int f = faces[n];
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        double area = mesh_->face_area(f);
       
        tmp = vf[n][i].Value(xf) - grad * (xf - xc);
        a1 += tmp * area;
        a2 += area;
      }

      uc[i](0, 0) = a1 / a2;
    } else if (order_ >= 2) {
      integrals_.poly().GetPolynomialCoefficients(v4);
      v4.Reshape(nd);
      uc[i](0, 0) = vdof(row) - (v4 * v5) / volume;
    }

    // calculate L2 projector
    if (type == Type::L2 && ndof_c > 0) {
      v5(0) = uc[i](0, 0);

      DG_Modal dg(order_, mesh_);
      dg.set_basis(TAYLOR_BASIS_NATURAL);

      DenseMatrix M, M2;
      DenseVector v6(nd - ndof_c);
      dg.MassMatrix(c, T, integrals_, M);

      M2 = M.SubMatrix(ndof_c, nd, 0, nd);
      M2.Multiply(v5, v6, false);

      for (int n = 0; n < ndof_c; ++n) {
        v4(n) = (*moments)(n) * mesh_->cell_volume(c);
      }

      for (int n = 0; n < nd - ndof_c; ++n) {
        v4(ndof_c + n) = v6(n);
      }

      M.Inverse();
      M.Multiply(v4, v5, false);

      uc[i].SetPolynomialCoefficients(v5);
      numi.ChangeBasisNaturalToRegular(c, uc[i]);
    }

    // set origin to zero
    AmanziGeometry::Point zero(d_);
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

