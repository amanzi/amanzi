/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MFD3D_ELECTROMAGNETICS_HH_
#define AMANZI_MFD3D_ELECTROMAGNETICS_HH_

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

class MFD3D_Electromagnetics : public MFD3D { 
 public:
  explicit MFD3D_Electromagnetics(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MFD3D(mesh) {};
  ~MFD3D_Electromagnetics() {};

  // required implementation of two consistency conditions
  // the inner product in the spave of edge-based functions is weighted by inverse(T)
  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);

  // consistency condition for stiffness matrix.
  // the inner product in the spave of edge-based functions is weighted by T.
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);

  int MassMatrix(int c, const Tensor& T, DenseMatrix& M);
  int MassMatrixOptimized(int c, const Tensor& T, DenseMatrix& W);

  int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W);
  int MassMatrixInverseOptimized(int c, const Tensor& T, DenseMatrix& W);

  int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A);
  int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C);
  int StiffnessMatrixExperimental(int c, const Tensor& T, DenseMatrix& A);

  // boundary and surface methods
  int L2consistencyBoundary(int f, const Tensor& K, DenseMatrix& R, DenseMatrix& Mf);
  int MassMatrixBoundary(int f, const Tensor& K, DenseMatrix& M);

 private:
  int L2consistency2D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  int L2consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);

  int L2consistencyInverse2D_(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc);
  int L2consistencyInverse3D_(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc);

  int H1consistency2DExperimental_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int H1consistency3DExperimental_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

