/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  High-order 3D Nedelec serendipity element of type 2: degrees of 
  freedom are moments on edges, faces and incide cell.
*/

#ifndef AMANZI_VEM_NEDELEC_SERENDIPITY_HH_
#define AMANZI_VEM_NEDELEC_SERENDIPITY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "GrammMatrix.hh"
#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "MFD3D.hh"
#include "Monomial.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class VEM_NedelecSerendipity : public DeRham_Edge { 
 public:
  VEM_NedelecSerendipity(const Teuchos::ParameterList& plist,
                              const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);
  ~VEM_NedelecSerendipity() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrix
  int L2consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  int MassMatrixFace(int f, const Tensor& K, DenseMatrix& M);

  // -- stiffness matrix
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // other methods
  int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C);
  void CurlMatrix(int c, DenseMatrix& C);

  // -- l2 projector
  void L2Cell(int c, const std::vector<VectorPolynomial>& ve,
              const std::vector<VectorPolynomial>& vf,
              const Polynomial* moments, VectorPolynomial& uc) {
    ProjectorCell_(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  void L2Face(int f, const std::vector<VectorPolynomial>& ve,
              const Polynomial* moments, VectorPolynomial& uf) {
    ProjectorFace_(f, ve, ProjectorType::L2, moments, uf);
  }

  // access
  PolynomialOnMesh& integrals() { return integrals_; }

  // support
  void CalculateDOFsOnBoundary(
      const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
      int c, const std::vector<VectorPolynomial>& ve,
      const std::vector<VectorPolynomial>& vf, DenseVector& vdof);

 protected:
  int L2consistency2D_(const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh,
                       int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& MG);

 private:
  void ProjectorCell_(const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
                      int c, const std::vector<VectorPolynomial>& ve,
                      const std::vector<VectorPolynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments, VectorPolynomial& uc);

  void ProjectorFace_(int f, const std::vector<VectorPolynomial>& ve,
                      const ProjectorType type,
                      const Polynomial* moments, VectorPolynomial& uf);

  // auxiliary functions
  void MatrixOfDofs_(
      int c, const Entity_ID_List& edges,
      const Basis_Regularized& basis,
      const NumericalIntegration& numi,
      DenseMatrix& N);

  void L2ProjectorsOnFaces_(
      int c, const Tensor& K, const Entity_ID_List& faces,
      std::vector<WhetStone::DenseMatrix>& vL2f, 
      std::vector<WhetStone::DenseMatrix>& vMGf,
      std::vector<Basis_Regularized>& vbasisf,
      std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> >& vcoordsys,
      int MGorder);

  void L2ProjectorOnEdge_(WhetStone::DenseMatrix& L2e, int order);

  void L2consistency3DFace_Method1_(
      const VectorPolynomial& p1,
      const VectorPolynomial& xyz,
      const SurfaceCoordinateSystem& coordsys,
      const Basis_Regularized& basis,
      const DenseMatrix& L2f,
      const DenseMatrix& MGf,
      DenseVector& p0v);

  void L2consistency3DFace_Method2_(
      int f,
      const VectorPolynomial& p1,
      const SurfaceCoordinateSystem& coordsys,
      const Basis_Regularized& basis,
      const DenseMatrix& L2f,
      const DenseMatrix& MGf,
      DenseVector& p0v);

 protected:
  using MFD3D::mesh_;
  using MFD3D::d_;

 private:
  int type_;
  PolynomialOnMesh integrals_;

 private:
  static RegisteredFactory<VEM_NedelecSerendipity> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
