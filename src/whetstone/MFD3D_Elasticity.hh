/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The package uses the formula M = Mc + Ms, where matrix Mc is build from a 
  consistency condition (Mc N = R) and matrix Ms is build from a stability 
  condition (Ms N = 0), to generate mass and stiffness matrices for a variety 
  of physics packages: flow, transport, thermal, and geomechanics. 
  The material properties are imbedded into the the matrix Mc. 

  Notation used below: M (mass), W (inverse of M), A (stiffness).
*/

#ifndef AMANZI_MFD3D_ELASTICITY_HH_
#define AMANZI_MFD3D_ELASTICITY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Elasticity : public virtual MFD3D { 
 public:
  MFD3D_Elasticity(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh),
      InnerProduct(mesh) {};
  ~MFD3D_Elasticity() {};

  // main method use edge-based DOFs (part of DeRham complex) 
  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) { return WHETSTONE_ELEMENTAL_MATRIX_OK; } 

  // -- inverse mass matrices
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) { return WHETSTONE_ELEMENTAL_MATRIX_OK; } 

  // -- steffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A);

  // overriding other methods
  int StiffnessMatrixOptimized(int c, const Tensor& T, DenseMatrix& A);
  int StiffnessMatrixMMatrix(int c, const Tensor& T, DenseMatrix& A);

  // complex DOFs
  // -- vectors at nodes, comal components on faces
  int H1consistencyBernardiRaugel(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);
  int StiffnessMatrixBernardiRaugel(int c, const Tensor& K, DenseMatrix& A);
  // -- divergence
  int DivergenceMatrixBernardiRaugel(int c, DenseMatrix& A);
  // -- advection
  int AdvectionMatrixBernardiRaugel(int c, DenseMatrix& A,
                                    const std::vector<AmanziGeometry::Point>& u);

 private:
  void MatrixMatrixProduct_(
      const DenseMatrix& A, const DenseMatrix& B, bool transposeB, DenseMatrix& AB);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

