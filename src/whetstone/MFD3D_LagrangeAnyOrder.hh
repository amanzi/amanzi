/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Lagrange-type element: degrees of freedom are nodal values and
  moments on edges, faces and inside cell.
*/

#ifndef AMANZI_MFD3D_LAGRANGE_ANY_ORDER_HH_
#define AMANZI_MFD3D_LAGRANGE_ANY_ORDER_HH_

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeAnyOrder : public MFD3D { 
 public:
  MFD3D_LagrangeAnyOrder(const Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_LagrangeAnyOrder() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override {
    Errors::Message msg("L2 consistency is not implemented for Lagrange_AnyOrder element.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override {
    Errors::Message msg("Mass matrix is not implemented for Lagrange_AnyOrder element.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override {
    if (d_ == 2) return H1consistency2D_(c, T, N, Ac);
    return H1consistency3D_(c, T, N, Ac);
  }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  virtual void H1Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, ProjectorType::H1, moments, uc);
  }

  virtual void H1Cell(int c, const DenseVector& dofs, Polynomial& uc) override {
    ProjectorCellFromDOFs_(c, dofs, ProjectorType::H1, uc);
  }

  // surface methods
  int H1consistencySurface(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int StiffnessMatrixSurface(int c, const Tensor& T, DenseMatrix& A);

  // access 
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  PolynomialOnMesh& integrals() { return integrals_; }

  // -- matrices that could be resused in other code
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  int H1consistency2D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int H1consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);

  void ProjectorCell_(int c, const std::vector<Polynomial>& vf,
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

}  // namespace WhetStone
}  // namespace Amanzi

#endif

