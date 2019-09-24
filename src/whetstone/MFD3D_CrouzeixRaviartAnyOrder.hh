/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#ifndef AMANZI_MFD3D_CROUZEIX_RAVIART_HO_HH_
#define AMANZI_MFD3D_CROUZEIX_RAVIART_HO_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MatrixPolynomial.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviartAnyOrder : public MFD3D { 
 public:
  MFD3D_CrouzeixRaviartAnyOrder(const Teuchos::ParameterList& plist,
                                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_CrouzeixRaviartAnyOrder() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override {
    Errors::Message msg("L2 consistency is not implemented for Crouzeix-Raviart space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override {
    Errors::Message msg("MassMatrix is not supported for Crouzeix-Raviart space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- projectors: base L2 and H1 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  virtual void H1Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, ProjectorType::H1, moments, uc);
  }

  // additional miscaleneous projectors
  void L2GradientCell(int c, const std::vector<VectorPolynomial>& vf,
                      const std::shared_ptr<DenseVector>& moments, MatrixPolynomial& uc) {
    ProjectorGradientCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  // generic code for multiple projectors
  void ProjectorCell_(int c, const std::vector<Polynomial>& vf,
                      const ProjectorType type, 
                      const Polynomial* moments, Polynomial& uc);

  void ProjectorGradientCell_(int c, const std::vector<VectorPolynomial>& vf,
                              const ProjectorType type, 
                              const std::shared_ptr<DenseVector>& moments, MatrixPolynomial& uc);

  // supporting routines
  void CalculateFaceDOFs_(int f, const Polynomial& vf, const Polynomial& pf,
                          DenseVector& vdof, int& row);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_CrouzeixRaviartAnyOrder> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

