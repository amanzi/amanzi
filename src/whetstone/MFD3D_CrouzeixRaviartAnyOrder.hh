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

#include "MeshLight.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "GrammMatrix.hh"
#include "MatrixObjects.hh"
#include "MFD3D.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviartAnyOrder : public MFD3D {
 public:
  MFD3D_CrouzeixRaviartAnyOrder(const Teuchos::ParameterList& plist,
                                const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  void L2GradientCell(int c,
                      const std::vector<VectorPolynomial>& vf,
                      const std::shared_ptr<DenseVector>& moments,
                      MatrixPolynomial& uc)
  {
    ProjectorGradientCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(mesh_, c, ve, vf, ProjectorType::H1, moments, uc);
  }

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  // generic code for multiple projectors
  void ProjectorCell_(const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh,
                      int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments,
                      Polynomial& uc);

  void ProjectorGradientCell_(int c,
                              const std::vector<VectorPolynomial>& vf,
                              const ProjectorType type,
                              const std::shared_ptr<DenseVector>& moments,
                              MatrixPolynomial& uc);

  // supporting routines
  void CalculateFaceDOFs_(int f,
                          const Polynomial& vf,
                          const Polynomial& pf,
                          DenseVector& vdof,
                          int& row);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_CrouzeixRaviartAnyOrder> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
