/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin modal method.
*/

#ifndef AMANZI_WHETSTONE_DG_MODAL_HH_
#define AMANZI_WHETSTONE_DG_MODAL_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

// Gauss quadrature on interval (0,1)
const double q1d_weights[4][4] = {
    1.0, 0.0, 0.0, 0.0,
    0.5, 0.5, 0.0, 0.0,
    0.277777777777778, 0.444444444444444, 0.277777777777778, 0.0,
    0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
};
const double q1d_points[4][4] = {
    0.5, 0.0, 0.0, 0.0,
    0.211324865405187, 0.788675134594813, 0.0, 0.0,
    0.112701665379258, 0.5, 0.887298334620742, 0.0,
    0.0694318442029737, 0.330009478207572, 0.669990521792428, 0.930568155797026
};

class DG_Modal { 
 public:
  DG_Modal(Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : order_(-1),
      mesh_(mesh),
      mesh0_(Teuchos::null),
      d_(mesh_->space_dimension()),
      method_(WHETSTONE_METHOD_VEM) {};

  DG_Modal(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
           Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : order_(-1),
      mesh_(mesh),
      mesh0_(mesh0),
      d_(mesh_->space_dimension()),
      method_(WHETSTONE_METHOD_VEM) {};

  DG_Modal(int order, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : order_(order), 
      mesh_(mesh),
      mesh0_(Teuchos::null),
      d_(mesh_->space_dimension()),
      method_(WHETSTONE_METHOD_VEM) {};

  ~DG_Modal() {};

  int MassMatrix(int c, const Tensor& K, DenseMatrix& M);
  int MassMatrix(int c, Polynomial& K, DenseMatrix& A);
  int AdvectionMatrixFace(int f, Polynomial& un, DenseMatrix& A);

  // support of advection schemes
  void FaceVelocity(int c, int f, std::vector<Polynomial>& v) const;
  Tensor FaceJacobian(int c, int f, const std::vector<Polynomial>& v,
                      const AmanziGeometry::Point& x) const;

  // finite elements use three meshes: reference, initial and target
  AmanziGeometry::Point FEM_Map(int c, const AmanziGeometry::Point& xref) const;
  Tensor FEM_Jacobian(int c, const AmanziGeometry::Point& xref) const;

  // miscalleneous
  void set_order(int order) { order_ = order; }
  void set_method(int method) { method_ = method; }

  // polynomial approximation of map x2 = F(x1)
  int LeastSquareFit(const std::vector<AmanziGeometry::Point>& x1, 
                     const std::vector<AmanziGeometry::Point>& x2,
                     std::vector<AmanziGeometry::Point>& u) const;

 private:
  void IntegrateMonomialsCell_(int c, Monomial& monomials);
  void IntegrateMonomialsFace_(int f, double factor, Monomial& monomials);
  void IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      double factor, Monomial& monomials);
  double IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      int ix, int iy, int jx, int jy,
      const std::vector<double>& factors, 
      const AmanziGeometry::Point& xc1, const AmanziGeometry::Point& xc2);

  // finite element maps 
  Tensor FEM_JacobianInternal_(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                               int c, const AmanziGeometry::Point& xref) const;
  Tensor FEM_FaceJacobian_(int c, int f, const std::vector<Polynomial>& v,
                           const AmanziGeometry::Point& xref) const;

  // virtual element maps 
  void VEM_FaceVelocity_(int c, int f, std::vector<Polynomial>& v) const;
  Tensor VEM_FaceJacobian_(int c, int f, const std::vector<Polynomial>& v,
                           const AmanziGeometry::Point& x) const;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;  // target mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;  // initial mesh 
  int order_, d_;
  int method_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

