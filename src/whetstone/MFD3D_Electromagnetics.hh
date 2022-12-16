/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The package uses the formula M = Mc + Ms, where matrix Mc is build from a
  consistency condition (Mc N = R) and matrix Ms is build from a stability
  condition (Ms N = 0), to generate mass and stiffness matrices for a variety
  of physics packages: flow, transport, thermal, and geomechanics.
  The material properties are imbedded into the the matrix Mc.

  Notation used below: M (mass), W (inverse of M), A (stiffness).
*/

#ifndef AMANZI_MFD3D_ELECTROMAGNETICS_HH_
#define AMANZI_MFD3D_ELECTROMAGNETICS_HH_


#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Electromagnetics : public DeRham_Edge {
 public:
  MFD3D_Electromagnetics(const Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
    : DeRham_Edge(mesh){};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override
  {
    if (d_ == 2)
      return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::NODE, DOF_Type::SCALAR, 1));
    else
      return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::EDGE, DOF_Type::SCALAR, 1));
  }

  // -- mass matrices
  // using InnerProductL2::MassMatrix;
  using DeRham_Edge::MassMatrix;
  using DeRham_Edge::MassMatrixInverse;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // other methods
  int MassMatrixOptimized(int c, const Tensor& T, DenseMatrix& M);
  int MassMatrixInverseOptimized(int c, const Tensor& T, DenseMatrix& M);

  int MassMatrixDiagonal(int c, const Tensor& T, DenseMatrix& M);

  using BilinearForm::StiffnessMatrix;
  int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C);
  int StiffnessMatrix_GradCorrection(int c, const Tensor& T, DenseMatrix& A);

  // curl matrix
  void CurlMatrix(int c, DenseMatrix& C);

  // boundary and surface methods
  int L2consistencyBoundary(int f, const Tensor& K, DenseMatrix& R, DenseMatrix& Mf);
  int MassMatrixBoundary(int f, const Tensor& K, DenseMatrix& M);

 private:
  int H1consistency2D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int H1consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);

 private:
  static RegisteredFactory<MFD3D_Electromagnetics> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
