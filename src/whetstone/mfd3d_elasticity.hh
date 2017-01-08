/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MFD3D_ELASTICITY_HH_
#define AMANZI_MFD3D_ELASTICITY_HH_

/*
  The package uses the formula M = Mc + Ms, where matrix Mc is build from a 
  consistency condition (Mc N = R) and matrix Ms is build from a stability 
  condition (Ms N = 0), to generate mass and stiffness matrices for a variety 
  of physics packages: flow, transport, thermal, and geomechanics. 
  The material properties are imbedded into the the matrix Mc. 

  Notation used below: M (mass), W (inverse of M), A (stiffness).
*/

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "mfd3d.hh"


namespace Amanzi {
namespace WhetStone {

class MFD3D_Elasticity : public MFD3D { 
 public:
  explicit MFD3D_Elasticity(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MFD3D(mesh) {};
  ~MFD3D_Elasticity() {};

  // Edges DOFs
  // -- consistency conditions
  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);

  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);

  int MassMatrix(int c, const Tensor& T, DenseMatrix& M) { return WHETSTONE_ELEMENTAL_MATRIX_OK; } 
  int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) { return WHETSTONE_ELEMENTAL_MATRIX_OK; } 

  int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A);
  int StiffnessMatrixOptimized(int c, const Tensor& T, DenseMatrix& A);
  int StiffnessMatrixMMatrix(int c, const Tensor& T, DenseMatrix& A);

  // complex sets of DOFs
  // -- vectors at nodes, comal components on faces
  int H1consistencyNode2Face1(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);

  int StiffnessMatrixNode2Face1(int c, const Tensor& K, DenseMatrix& A);

 private:
  void MatrixMatrixProduct_(
      const DenseMatrix& A, const DenseMatrix& B, bool transposeB, DenseMatrix& AB);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

