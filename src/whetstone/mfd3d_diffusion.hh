/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_WHETSTONE_MFD3D_DIFFUSION_HH_
#define AMANZI_WHETSTONE_MFD3D_DIFFUSION_HH_

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

class MFD3D_Diffusion : public MFD3D { 
 public:
  explicit MFD3D_Diffusion(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MFD3D(mesh) {};
  ~MFD3D_Diffusion() {};

  // basic mimetic discretization methods use permeability tensor K
  // the inner product in the spave of face-based functions is weighted by
  // inverse of K. 
  int L2consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  int L2consistencyInverse(int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);

  int H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);
  int H1consistencyEdge(int cell, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);

  int MassMatrix(int c, const Tensor& K, DenseMatrix& M);
  int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseOptimized(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseMMatrixHex(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseMMatrix(int c, const Tensor& K, DenseMatrix& W);

  int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A);
  int StiffnessMatrixOptimized(int c, const Tensor& K, DenseMatrix& A);
  int StiffnessMatrixMMatrix(int c, const Tensor& K, DenseMatrix& A);

  int StiffnessMatrixEdge(int c, const Tensor& K, DenseMatrix& A);

  // natural scaling of fluxes which was found to be the right way.
  int L2consistencyInverseScaled(int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc);
  int MassMatrixInverseScaled(int c, const Tensor& K, DenseMatrix& W); 
  int MassMatrixInverseOptimizedScaled(int c, const Tensor& K, DenseMatrix& W);

  // experimental methods
  // -- tensor is product k K
  int L2consistencyInverseDivKScaled(int c, const Tensor& K,
                                     double kmean, const AmanziGeometry::Point& kgrad,
                                     DenseMatrix& R, DenseMatrix& Wc);
  int MassMatrixInverseDivKScaled(int c, const Tensor& K,
                                  double kmean, const AmanziGeometry::Point& kgrad, DenseMatrix& W);

  // -- non-symmetric tensor K
  int MassMatrixNonSymmetric(int c, const Tensor& K, DenseMatrix& M);
  int MassMatrixInverseNonSymmetric(int c, const Tensor& K, DenseMatrix& W);

  // surface methods
  int L2consistencyInverseSurface(int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc);
  int MassMatrixInverseSurface(int c, const Tensor& K, DenseMatrix& W);

  // primary related discetization methods
  int MassMatrixInverseSO(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseTPFA(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseDiagonal(int c, const Tensor& K, DenseMatrix& W);

  // a posteriori error estimate
  int RecoverGradient_MassMatrix(int c, const std::vector<double>& solution, 
                                 AmanziGeometry::Point& gradient);

  int RecoverGradient_StiffnessMatrix(int c, const std::vector<double>& solution, 
                                      AmanziGeometry::Point& gradient);

  // utils
  double Transmissibility(int f, int c, const Tensor& K);

 private:  
  // stability methods (add stability matrix, M += Mstab)
  int StabilityMMatrixHex_(int c, const Tensor& K, DenseMatrix& M);
  void RescaleMassMatrixInverse_(int c, DenseMatrix& W);
  void StabilityScalarNonSymmetric_(int c, DenseMatrix& N, DenseMatrix& M);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

