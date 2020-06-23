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
    (3) moments of edges, groupped by edge (in 3D);
    (4) moments inside cell.
*/

#ifndef AMANZI_MFD3D_LAGRANGE_ANY_ORDER_HH_
#define AMANZI_MFD3D_LAGRANGE_ANY_ORDER_HH_

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "GrammMatrix.hh"
#include "MFD3D.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeAnyOrder : public MFD3D { 
 public:
  MFD3D_LagrangeAnyOrder(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh) {};
  MFD3D_LagrangeAnyOrder(const Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_LagrangeAnyOrder() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override {
    Errors::Message msg("L2 consistency is not implemented for LagrangeAnyOrder element.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override {
    Errors::Message msg("Mass matrix is not implemented for LagrangeAnyOrder element.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override {
    if (d_ == 2) return H1consistency2D_<AmanziMesh::Mesh>(mesh_, c, T, N, Ac);
    return H1consistency3D_(c, T, N, Ac, true);
  }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, ve, vf, ProjectorType::L2, moments, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, ve, vf, ProjectorType::H1, moments, uc);
  }

  virtual void H1Cell(int c, const DenseVector& dofs, Polynomial& uc) override {
    ProjectorCellFromDOFs_(c, dofs, ProjectorType::H1, uc);
  }

  // surface methods
  int StiffnessMatrixSurface(int c, const Tensor& T, DenseMatrix& A);

  // access 
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  PolynomialOnMesh& integrals() { return integrals_; }

  // -- matrices that could be resused in other code
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 protected:
  template<typename MyMesh>
  int H1consistency2D_(const Teuchos::RCP<const MyMesh>& mymesh,
                       int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);

  int H1consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac, bool doAc);

 private:
  void ProjectorCell_(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments, Polynomial& uc);

  void ProjectorCellFromDOFs_(int c, const DenseVector& dofs,
                              const ProjectorType type, Polynomial& uc);

  std::vector<Polynomial> ConvertMomentsToPolynomials_(int order);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_LagrangeAnyOrder> factory_;
};


/* ******************************************************************
* High-order consistency condition for the stiffness matrix. 
****************************************************************** */
template <class MyMesh>
int MFD3D_LagrangeAnyOrder::H1consistency2D_(
    const Teuchos::RCP<const MyMesh>& mymesh,
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mymesh->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mymesh->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mymesh->cell_centroid(c); 
  double volume = mymesh->cell_volume(c); 

  // calculate degrees of freedom 
  Polynomial poly(d, order_), pf, pc;
  if (order_ > 1) {
    pf.Reshape(d - 1, order_ - 2);
    pc.Reshape(d, order_ - 2);
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
  NumericalIntegration<MyMesh> numi(mymesh);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_ - 2, integrals_);

  // selecting regularized basis
  Basis_Regularized<MyMesh> basis;
  basis.Init(mymesh, c, order_, integrals_.poly());

  // populate matrices N and R
  R_.PutScalar(0.0);
  N.PutScalar(0.0);

  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) { 
    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d, index, factor);
    cmono.set_origin(xc);  

    // N: degrees of freedom at vertices
    auto grad = Gradient(cmono);
     
    polys[0] = &cmono;

    int col = it.PolynomialPosition();
    int row(nnodes);

    AmanziGeometry::Point xv(d);
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      mymesh->node_get_coordinates(v, &xv);
      N(i, col) = cmono.Value(xv);
    }

    // N and R: degrees of freedom on faces 
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      double area = mymesh->face_area(f);
      const AmanziGeometry::Point& xf = mymesh->face_centroid(f); 
      AmanziGeometry::Point normal = mymesh->face_normal(f);

      // local coordinate system with origin at face centroid
      auto coordsys = std::make_shared<SurfaceCoordinateSystem>(xf, normal);

      normal *= dirs[i];
      AmanziGeometry::Point conormal = K * normal;

      Entity_ID_List face_nodes;
      mymesh->face_get_nodes(f, &face_nodes);
      int nfnodes = face_nodes.size();

      if (order_ == 1 && col > 0) {
        for (int j = 0; j < nfnodes; j++) {
          int v = face_nodes[j];
          int pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
          R_(pos, col) += factor * conormal[col - 1] / 2;
        }
      } else if (col > 0) {
        int v, pos0, pos1;
        AmanziGeometry::Point x0(d), x1(d), xm(d), sm(d);

        Polynomial tmp = grad * conormal;

        v = face_nodes[0];
        pos0 = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
        mymesh->node_get_coordinates(v, &x0);

        v = face_nodes[1];
        pos1 = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
        mymesh->node_get_coordinates(v, &x1);

        if (order_ == 2) {
          // Simpson rule with 3 points
          double q0 = tmp.Value(x0);
          double q1 = tmp.Value(x1);
          double qmid = tmp.Value(mymesh->face_centroid(f));

          R_(pos0, col) += (q0 - qmid) / 6;
          R_(pos1, col) += (q1 - qmid) / 6;
          R_(row,  col) = qmid;
        } else if (order_ > 2) {
          if (col < 3) {
            // constant gradient contributes only to 0th moment 
            R_(row, col) += tmp(0);
          } else {
            auto polys_f = ConvertMomentsToPolynomials_(order_);

            // Gauss-Legendre quadrature rule with (order_) points
            int m(order_ - 1); 
            for (int n = 0; n < order_; ++n) { 
              xm = x0 * q1d_points[m][n] + x1 * (1.0 - q1d_points[m][n]);
              sm[0] = 0.5 - q1d_points[m][n];

              factor = q1d_weights[m][n] * tmp.Value(xm);
              R_(pos0, col) += polys_f[0].Value(sm) * factor;
              R_(pos1, col) += polys_f[1].Value(sm) * factor;

              for (int k = 0; k < m; ++k) { 
                R_(row + k, col) += polys_f[k + 2].Value(sm) * factor;
              }
            }
          }
        }
      }

      if (order_ > 1) {
        for (auto jt = pf.begin(); jt < pf.end(); ++jt) {
          const int* jndex = jt.multi_index();
          Polynomial fmono(d - 1, jndex, 1.0);
          fmono.InverseChangeCoordinates(xf, *coordsys->tau());  

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
        for (int i = 0; i < d; ++i) {
          multi_index[i] = index[i] + jndex[i];
          nm += multi_index[i];
        }

        int m = MonomialSetPosition(d, multi_index);
        factor = basis.monomial_scales()[it.MonomialSetOrder()] *
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

}  // namespace WhetStone
}  // namespace Amanzi

#endif

