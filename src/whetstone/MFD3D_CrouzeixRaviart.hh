/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#ifndef AMANZI_MFD3D_CROUZEIX_RAVIART_HH_
#define AMANZI_MFD3D_CROUZEIX_RAVIART_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviart : public virtual MFD3D { 
 public:
  MFD3D_CrouzeixRaviart(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh),
      InnerProduct(mesh),
      order_(1) {};
  ~MFD3D_CrouzeixRaviart() {};

  // required methods
  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) { return -1; }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) { return -1; } 

  // -- inverse mass matrices
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) { return -1; }
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& M) { return -1; } 

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) {;
    if (order_ == 1) {
      return H1consistencyLO_(c, T, N, Ac);
    } else {
      DenseMatrix R, G;
      return H1consistencyHO(c, order_, T, N, R, G, Ac);
    }
  }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) {
    if (order_ == 1) {
      return StiffnessMatrixLO_(c, T, A);
    } else {
      DenseMatrix R, G;
      return StiffnessMatrixHO(c, order_, T, R, G, A);
    }
  }

  // -- projectors
  virtual void H1Cell(
      int c, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) override {
    ProjectorCell_HO_(c, vf, Type::H1, false, moments, uc);
  }

  // harmonic projector calculates and returns cell-moments
  void H1CellHarmonic(
      int c, const std::vector<VectorPolynomial>& vf,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) {
    if (order_ == 1) {
      ProjectorCell_LO_(c, vf, uc);
    } else {
      ProjectorCell_HO_(c, vf, Type::H1, true, moments, uc);
    }
  }

  void H1FaceHarmonic(
      int f, const AmanziGeometry::Point& p0,
      const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const;

  // high-order methods
  int H1consistencyHO(int c, int order, const Tensor& T,
                      DenseMatrix& N, DenseMatrix& R, DenseMatrix& G, DenseMatrix& Ac);
  int StiffnessMatrixHO(int c, int order, const Tensor& T,
                        DenseMatrix& R, DenseMatrix& G, DenseMatrix& A);

  // miscalleneous
  void set_order(int order) { order_ = order; }

  // access 
  // -- integrals of monomials in high-order schemes could be reused
  const Polynomial& integrals() const { return integrals_; }

 private:
  // efficient implementation of low-order methods
  int H1consistencyLO_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int StiffnessMatrixLO_(int c, const Tensor& T, DenseMatrix& A);

  // simple implementation of low-order elliptic projectors
  void ProjectorCell_LO_(
      int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc);

  void ProjectorCell_HO_(
      int c, const std::vector<VectorPolynomial>& vf,
      const Projectors::Type type, bool is_harmonic,
      const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc);

 private:
  int order_;
  Polynomial integrals_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

