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
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class DG_Modal : public BilinearForm {
 public:
  DG_Modal(Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : numi_(0, mesh),
      mesh_(mesh),
      order_(-1),
      basis_(TAYLOR_BASIS_NORMALIZED),
      d_(mesh_->space_dimension()) {};

  DG_Modal(int order, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : numi_(order, mesh),
      order_(order), 
      mesh_(mesh),
      basis_(TAYLOR_BASIS_NORMALIZED),
      d_(mesh_->space_dimension()) {};

  ~DG_Modal() {};

  // required member functions
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
  
  // create polynomial given array of its coefficients
  Polynomial CalculatePolynomial(int c, const std::vector<double>& coefs) const;

  // placeholder of elliptic projector
  void CoVelocityCell(int c, const std::vector<const Polynomial*> vf, VectorPolynomial& vc);

  // miscalleneous
  void set_order(int order) { order_ = order; }
  void set_basis(int basis) { basis_ = basis; }

 private:
  void UpdateIntegrals_(int c, int order);
  void UpdateScales_(int c, int order);
  void ChangeBasis_(int c, DenseMatrix& A);
  void ChangeBasis_(int c1, int c2, DenseMatrix& A);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  NumericalIntegration numi_;
  int order_, d_;
  int basis_;

  VectorPolynomial integrals_;  // integrals of non-normalized monomials
  VectorPolynomial scales_a_;   // partial orthonormalization of Taylor basis
  VectorPolynomial scales_b_;  
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

