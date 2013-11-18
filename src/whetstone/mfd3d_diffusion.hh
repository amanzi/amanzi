/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Version: 2.0
  Release name: naka-to.
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
#include "tensor.hh"
#include "mfd3d.hh"


namespace Amanzi {
namespace WhetStone {

class MFD3D_Diffusion : public MFD3D { 
 public:
  explicit MFD3D_Diffusion(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MFD3D(mesh) {};
  ~MFD3D_Diffusion() {};

  int L2consistency(int cell, const Tensor& T,
                    DenseMatrix& N, DenseMatrix& Mc);

  int L2consistencyInverse(int cell, const Tensor& permeability,
                           DenseMatrix& R, DenseMatrix& Wc);

  int H1consistency(int cell, const Tensor& T, 
                    DenseMatrix& N, DenseMatrix& Mc);

  int MassMatrix(int cell, const Tensor& permeability, DenseMatrix& M);

  int MassMatrixInverse(int cell, const Tensor& permeability, DenseMatrix& W);

  int StiffnessMatrix(int cell, const Tensor& permeability, DenseMatrix& A);

  int StiffnessMatrixMMatrix(int cell, const Tensor& permeability, DenseMatrix& A);

  // experimental methods
  int L2consistencyInverseScaled(int cell, const Tensor& permeability,
                                 DenseMatrix& R, DenseMatrix& Wc);

  int MassMatrixInverseScaled(int cell, const Tensor& permeability,
                              DenseMatrix& W);

  int MassMatrixInverseOptimized(int cell, const Tensor& permeability,
                                 DenseMatrix& W);

  int MassMatrixInverseOptimizedScaled(int cell, const Tensor& permeability,
                                       DenseMatrix& W);

  // primary related discetization methods
  int MassMatrixInverseMMatrixHex(int cell, const Tensor& permeability, DenseMatrix& W);
  int MassMatrixInverseMMatrix(int cell, const Tensor& permeability, DenseMatrix& W);

  int MassMatrixInverseSO(int cell, const Tensor& permeability, DenseMatrix& W);

  int MassMatrixInverseTPFA(int cell, const Tensor& permeability, DenseMatrix& W);

  int MassMatrixInverseDiagonal(int cell, const Tensor& permeability, DenseMatrix& W);

  // a posteriori error estimate
  int RecoverGradient_MassMatrix(int cell,
                                 const std::vector<double>& solution, 
                                 AmanziGeometry::Point& gradient);

  int RecoverGradient_StiffnessMatrix(int cell,
                                      const std::vector<double>& solution, 
                                      AmanziGeometry::Point& gradient);

  // access
  double simplex_functional() { return simplex_functional_; }
  int simplex_num_itrs() { return simplex_num_itrs_; }

 private:  
  // stability methods (add matrix Ms in M = Mc + Ms)
  int StabilityMMatrixHex_(int cell, const Tensor& T, DenseMatrix& Mc, DenseMatrix& M);

  int StabilityMMatrix_(int cell, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& M);

  void RescaleMassMatrixInverse_(int cell, DenseMatrix& W);

  int SimplexFindFeasibleSolution_(DenseMatrix& T, int m1, int m2, int m3, int* izrow, int* iypos);
  void SimplexPivotElement_(DenseMatrix& T, int kp, int* ip);
  void SimplexExchangeVariables_(DenseMatrix& T, int kp, int ip);

  double simplex_functional_;
  int simplex_num_itrs_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

