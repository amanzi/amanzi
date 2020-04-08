/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  High-order 3D Raviart-Thomas serendipity element.
*/

#ifndef AMANZI_VEM_RAVIART_THOMAS_SERENDIPITY_HH_
#define AMANZI_VEM_RAVIART_THOMAS_SERENDIPITY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "GrammMatrix.hh"
#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "FunctionPower.hh"
#include "MFD3D.hh"
#include "Monomial.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "SurfaceMiniMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class VEM_RaviartThomasSerendipity : public MFD3D {
 public:
  VEM_RaviartThomasSerendipity(const Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~VEM_RaviartThomasSerendipity() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrix
  virtual int L2consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override;
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override { return 0; }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override { return 0; }

  // -- l2 projector
  void L2Cell(int c, const std::vector<VectorPolynomial>& ve,
              const std::vector<VectorPolynomial>& vf,
              const Polynomial* moments, VectorPolynomial& uc) {};

  // access
  DenseMatrix& G() { return G_; }
  const std::vector<DenseMatrix>& MGf() { return MGf_; }

 protected:
  using MFD3D::mesh_;
  using MFD3D::d_;

 private:
  PolynomialOnMesh integrals_;
  DenseMatrix G_;

  bool save_face_matrices_;
  std::vector<DenseMatrix> MGf_;

 private:
  static RegisteredFactory<VEM_RaviartThomasSerendipity> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
