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

#ifndef AMANZI_MFD3D_ELECTROMAGNETICS_HH_
#define AMANZI_MFD3D_ELECTROMAGNETICS_HH_


#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Electromagnetics : public virtual MFD3D,
                               public virtual DeRham_Edge {
 public:
  MFD3D_Electromagnetics(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh),
      InnerProduct(mesh) {};
  ~MFD3D_Electromagnetics() {};

  // main methods
  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A);

  // -- other matrices
  virtual int DivergenceMatrix(int c, DenseMatrix& A) { return -1; }

  // other methods
  int MassMatrixOptimized(int c, const Tensor& T, DenseMatrix& M);
  int MassMatrixInverseOptimized(int c, const Tensor& T, DenseMatrix& M);

  int MassMatrixDiagonal(int c, const Tensor& T, DenseMatrix& M);

  int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C);
  int StiffnessMatrixExperimental(int c, const Tensor& T, DenseMatrix& A);

  // boundary and surface methods
  int L2consistencyBoundary(int f, const Tensor& K, DenseMatrix& R, DenseMatrix& Mf);
  int MassMatrixBoundary(int f, const Tensor& K, DenseMatrix& M);

 private:
  int H1consistency2DExperimental_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int H1consistency3DExperimental_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

