/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __MFD3D_ELASTICITY_HH__
#define __MFD3D_ELASTICITY_HH__

/*
This is the discretization package.

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

class MFD3D_Elasticity : public MFD3D { 
 public:
  explicit MFD3D_Elasticity(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : MFD3D(mesh) {};
  ~MFD3D_Elasticity() {};

  // required implementation of two consistency conditions
  int L2consistency(int cell, const Tensor& deformation,
                    Teuchos::SerialDenseMatrix<int, double>& N,
                    Teuchos::SerialDenseMatrix<int, double>& Mc);

  int L2consistencyInverse(int cell, const Tensor& deformation,
                           Teuchos::SerialDenseMatrix<int, double>& R,
                           Teuchos::SerialDenseMatrix<int, double>& Wc) { return WHETSTONE_ELEMENTAL_MATRIX_OK; }

  int H1consistency(int cell, const Tensor& deformation,
                    Teuchos::SerialDenseMatrix<int, double>& N,
                    Teuchos::SerialDenseMatrix<int, double>& Mc);

  int MassMatrix(int cell, const Tensor& deformation,
                 Teuchos::SerialDenseMatrix<int, double>& M) { return WHETSTONE_ELEMENTAL_MATRIX_OK; } 

  int MassMatrixInverse(int cell, const Tensor& deformation,
                        Teuchos::SerialDenseMatrix<int, double>& W) { return WHETSTONE_ELEMENTAL_MATRIX_OK; } 

  int StiffnessMatrix(int cell, const Tensor& deformation,
                      Teuchos::SerialDenseMatrix<int, double>& A);

 private:
  void MatrixMatrixProduct_(
      const Teuchos::SerialDenseMatrix<int, double>& A,
      const Teuchos::SerialDenseMatrix<int, double>& B, bool transposeB,
      Teuchos::SerialDenseMatrix<int, double>& AB);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

