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

#ifndef AMANZI_MFD3D_LAGRANGE_HH_
#define AMANZI_MFD3D_LAGRANGE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Lagrange : public MFD3D { 
 public:
  MFD3D_Lagrange(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh),
      InnerProduct(mesh) {};
  ~MFD3D_Lagrange() {};

  // required methods
  // -- mass matrices
  virtual int L2consistency(
      int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override { return -1; }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override { return -1; } 

  // -- inverse mass matrices
  virtual int L2consistencyInverse(
     int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) override { return -1; }
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& M) override { return -1; } 

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- projectors
  virtual void L2Cell(
      int c, const std::vector<Polynomial>& vf,
      Polynomial& moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, Type::L2, moments, uc);
  }

  virtual void H1Cell(
      int c, const std::vector<Polynomial>& vf,
      Polynomial& moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, Type::H1, moments, uc);
  }

  // access 
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  PolynomialOnMesh& integrals() { return integrals_; }

  // -- matrices that could be resused in other code
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  void ProjectorCell_(
      int c, const std::vector<Polynomial>& vf,
      const Projectors::Type type,
      Polynomial& moments, Polynomial& uc);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

