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

#include "BilinearForm.hh"
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

class DG_Modal : public BilinearForm { 
 public:
  DG_Modal(Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : order_(-1),
      basis_(TAYLOR_BASIS_NORMALIZED),
      mesh_(mesh),
      d_(mesh_->space_dimension()) {};

  DG_Modal(int order, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : order_(order), 
      basis_(TAYLOR_BASIS_NORMALIZED),
      mesh_(mesh),
      d_(mesh_->space_dimension()) {};

  ~DG_Modal() {};

  // requires member functions
  // -- mass matrices
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& M);
  virtual int MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M);

  // -- stiffness matrices (coming soon)
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) { return 0; }
  virtual int StiffnessMatrixPoly(int c, const Polynomial& K, DenseMatrix& A) { return 0; }

  // -- advection matrices
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A, bool grad_on_test) { return 0; }
  virtual int AdvectionMatrixPoly(int c, const VectorPolynomial& uc, DenseMatrix& A, bool grad_on_test);
  int FluxMatrixPoly(int f, const Polynomial& uf, DenseMatrix& A, bool jump_on_test);

  // interfaces that are not used
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) { return 0; }
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) { return 0; }
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc) { return 0; }

  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) { return 0; }
  virtual int DivergenceMatrix(int c, DenseMatrix& A) { return 0; }

  // scaling of Taylor basis function: \psi_k -> a (\psi_k - b \psi_0)
  void TaylorBasis(int c, const Iterator& it, double* a, double* b);
  
  // calculate polynomial given arrsy of its coefficients
  Polynomial CalculatePolynomial(int c, const std::vector<double>& coefs) const;

  // integration tools
  double IntegratePolynomialFace(int f, const Polynomial& poly) const {
    const std::vector<const Polynomial*> polys(1, &poly);
    return IntegratePolynomialsFace_(f, polys);
  }

  double IntegratePolynomialEdge(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const Polynomial& poly) const {
    const std::vector<const Polynomial*> polys(1, &poly);
    return IntegratePolynomialsEdge_(x1, x2, polys);
  }

  // miscalleneous
  void set_order(int order) { order_ = order; }
  void set_basis(int basis) { basis_ = basis; }

 private:
  // specialized routines optimized for non-normalized Taylor basis
  void IntegrateMonomialsCell_(int c, Monomial& monomials);
  void IntegrateMonomialsFace_(int c, int f, double factor, Monomial& monomials);
  void IntegrateMonomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      double factor, Monomial& monomials);

  // integraton of a product of polynomials with potentialy different origins
  double IntegratePolynomialsFace_(
      int f, const std::vector<const Polynomial*>& polys) const;

  double IntegratePolynomialsEdge_(
      const AmanziGeometry::Point& x1, const AmanziGeometry::Point& x2,
      const std::vector<const Polynomial*>& polys) const;

 private:
  void UpdateIntegrals_(int c, int order);
  void UpdateScales_(int c, int order);
  void ChangeBasis_(int c, DenseMatrix& A);
  void ChangeBasis_(int c1, int c2, DenseMatrix& A);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int order_, d_;
  int basis_;

  VectorPolynomial integrals_;  // integrals of non-normalized monomials
  VectorPolynomial scales_a_;   // partial orthonormalization of Taylor basis
  VectorPolynomial scales_b_;  
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

