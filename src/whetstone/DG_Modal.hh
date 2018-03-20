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

// Amanzi
#include "Mesh.hh"
#include "Point.hh"

// WhetStone
#include "BilinearForm.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class DG_Modal : public BilinearForm {
 public:
  DG_Modal(Teuchos::RCP<const AmanziMesh::Mesh> mesh) 
    : numi_(mesh),
      mesh_(mesh),
      order_(-1),
      basis_(TAYLOR_BASIS_NORMALIZED_ORTHO),
      d_(mesh_->space_dimension()) {};

  DG_Modal(int order, Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : numi_(mesh),
      order_(order), 
      mesh_(mesh),
      basis_(TAYLOR_BASIS_NORMALIZED_ORTHO),
      d_(mesh_->space_dimension()) {};

  ~DG_Modal() {};

  // basic member functions
  // -- mass matrices
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& M) override;
  virtual int MassMatrix(int c, const Polynomial& K, DenseMatrix& M) override;
  int MassMatrix(int c, const Tensor& K, PolynomialOnMesh& integrals, DenseMatrix& M);

  // -- inverse mass matrices
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W) override {
    int ok = MassMatrix(c, K, W);
    W.Inverse();
    return ok;
  }

  virtual int MassMatrixInverse(int c, const Polynomial& K, DenseMatrix& W) override {
    int ok = MassMatrix(c, K, W);
    W.Inverse();
    return ok;
  }

  // -- stiffness matrices
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // -- advection matrices
  virtual int AdvectionMatrix(int c, const VectorPolynomial& uc, DenseMatrix& A, bool grad_on_test) override;

  // -- flux matrices
  int FluxMatrixUpwind(int f, const Polynomial& uf, DenseMatrix& A, bool jump_on_test);
  int FluxMatrixRusanov(int f, const VectorPolynomial& uc1, const VectorPolynomial& uc2,
                        const Polynomial& uf, DenseMatrix& A);

  // -- interface matrices: jumps and penalty
  int FaceMatrixJump(int f, const Tensor& K1, const Tensor& K2, DenseMatrix& A);
  int FaceMatrixPenalty(int f, double Kf, DenseMatrix& A);

  // interfaces that are not used
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override { return 0; }
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) override { return 0; }
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc) override { return 0; }

  // scaling of Taylor basis function: \psi_k -> a (\psi_k - b \psi_0)
  void TaylorBasis(int c, const Iterator& it, double* a, double* b);
  
  // create polynomial given array of its coefficients
  Polynomial CalculatePolynomial(int c, const DenseVector& coefs) const;

  // miscalleneous
  // -- order of polynomials in each cell
  void set_order(int order) { order_ = order; }
  // -- modifications of Taylor basis. Available options are natural, 
  //    normalized by volume, normalized aind orthogonalized to constant
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

