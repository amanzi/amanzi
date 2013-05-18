/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: aka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __MFD3D_DIFFUSION_HH__
#define __MFD3D_DIFFUSION_HH__

/*
This is the discretization package, release alpha.

The package uses the formula M = Mc + Ms, where matrix Mc is build from a 
consistency condition (Mc N = R) and matrix Ms is build from a stability 
condition (Ms N = 0), to generate mass and stiffness matrices for a variety 
of physics packages: flow, transport, thermal, and geomechanics. 
The material properties are imbedded into the the matrix Mc. 

Notation used below: M (mass), W (inverse of M), A (stiffness).

IMPORTANT: all matrices must be reshaped before calling member functions.
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "tensor.hh"
#include "mfd3d.hh"


namespace Amanzi {
namespace WhetStone {

class MFD3D_Diffusion : public MFD3D { 
 public:
  explicit MFD3D_Diffusion(Teuchos::RCP<const AmanziMesh::Mesh> mesh) :MFD3D(mesh) {};
  ~MFD3D_Diffusion() {};

  int L2consistency(int cell, const Tensor& T,
                    Teuchos::SerialDenseMatrix<int, double>& N,
                    Teuchos::SerialDenseMatrix<int, double>& Mc);

  int L2consistencyInverse(int cell, const Tensor& permeability,
                           Teuchos::SerialDenseMatrix<int, double>& R,
                           Teuchos::SerialDenseMatrix<int, double>& Wc);

  int H1consistency(int cell, const Tensor& T,
                    Teuchos::SerialDenseMatrix<int, double>& N,
                    Teuchos::SerialDenseMatrix<int, double>& Mc);

  int MassMatrix(int cell, const Tensor& permeability,
                 Teuchos::SerialDenseMatrix<int, double>& M);

  int MassMatrixInverse(int cell, const Tensor& permeability,
                        Teuchos::SerialDenseMatrix<int, double>& W);

  int StiffnessMatrix(int cell, const Tensor& permeability,
                      Teuchos::SerialDenseMatrix<int, double>& A);

  // experimental methods
  int L2consistencyInverseScaled(int cell, const Tensor& permeability,
                                 Teuchos::SerialDenseMatrix<int, double>& R,
                                 Teuchos::SerialDenseMatrix<int, double>& Wc);

  int MassMatrixInverseScaled(int cell, const Tensor& permeability,
                              Teuchos::SerialDenseMatrix<int, double>& W);

  int MassMatrixInverseOptimized(int cell, const Tensor& permeability,
                                 Teuchos::SerialDenseMatrix<int, double>& W);

  int MassMatrixInverseOptimizedScaled(int cell, const Tensor& permeability,
                                       Teuchos::SerialDenseMatrix<int, double>& W);

  // primary related discetization methods
  int MassMatrixInverseHex(int cell, const Tensor& permeability,
                           Teuchos::SerialDenseMatrix<int, double>& W);

  int MassMatrixInverseSO(int cell, const Tensor& permeability,
                          Teuchos::SerialDenseMatrix<int, double>& W);

  int MassMatrixInverseDiagonal(int cell, const Tensor& permeability,
                                Teuchos::SerialDenseMatrix<int, double>& W);

  // a posteriori error estimate
  int RecoverGradient_MassMatrix(int cell,
                                 const std::vector<double>& solution, 
                                 AmanziGeometry::Point& gradient);

  int RecoverGradient_StiffnessMatrix(int cell,
                                      const std::vector<double>& solution, 
                                      AmanziGeometry::Point& gradient);


 private:  
  // supporting stability methods (add matrix Ms in M = Mc + Ms)
  int StabilityMonotoneHex(int cell, const Tensor& T,
                           Teuchos::SerialDenseMatrix<int, double>& Mc,
                           Teuchos::SerialDenseMatrix<int, double>& M);

 private:
  void RescaleMassMatrixInverse_(int cell, Teuchos::SerialDenseMatrix<int, double>& W);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

