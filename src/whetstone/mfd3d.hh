/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Version: 2.0
Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __MFD3D_HH__
#define __MFD3D_HH__

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

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "tensor.hh"

namespace Amanzi {
namespace WhetStone {

const int WHETSTONE_ELEMENTAL_MATRIX_OK = 0;
const int WHETSTONE_ELEMENTAL_MATRIX_WRONG = 1;
const int WHETSTONE_ELEMENTAL_MATRIX_PASSED = 2;
const int WHETSTONE_ELEMENTAL_MATRIX_FAILED = 4;  // only for unexpected situations

const int WHETSTONE_STABILITY_GENERIC = 1;
const int WHETSTONE_STABILITY_GENERIC_SCALED = 2;
const int WHETSTONE_STABILITY_OPTIMIZED_DMP = 3;
const int WHETSTONE_STABILITY_OPTIMIZED_GEOMETRY = 4;

const int WHETSTONE_MAX_SPATIAL_DIMENSION = 3;

class MFD3D { 
 public:
  explicit MFD3D(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D() {};

  virtual int L2consistency(int cell, const Tensor& T,
                            DenseMatrix& N, DenseMatrix& Mc) = 0;

  virtual int L2consistencyInverse(int cell, const Tensor& T,
                                   DenseMatrix& R, DenseMatrix& Wc) = 0;

  virtual int H1consistency(int cell, const Tensor& T,
                            DenseMatrix& N, DenseMatrix& Mc) = 0;

  virtual int MassMatrix(int cell, const Tensor& T, DenseMatrix& M) = 0; 

  virtual int MassMatrixInverse(int cell, const Tensor& T, DenseMatrix& W) = 0; 

  virtual int StiffnessMatrix(int cell, const Tensor& T, DenseMatrix& A) = 0; 

  // experimental methods (for stability region analysis; unit test)
  double CalculateStabilityScalar(DenseMatrix& Mc);
  void ModifyStabilityScalingFactor(double factor);

  // access members
  double scaling_factor() { return scaling_factor_; }
  double scalar_stability() { return scalar_stability_; }

  // expension of mesh API (must be removed from this class lipnikov@lanl.gov)
  int cell_get_face_adj_cell(const int cell, const int face);

 protected:
  // supporting stability methods (add matrix Ms in M = Mc + Ms)
  void StabilityScalar(int cell, DenseMatrix& N,  // use R, Wc, W for the inverse matrix
                       DenseMatrix& Mc, DenseMatrix& M);

  int StabilityOptimized(const Tensor& T, DenseMatrix& N, 
                         DenseMatrix& Mc, DenseMatrix& M);

  int StabilityMonotoneHex(int cell, const Tensor& T,
                           DenseMatrix& Mc, DenseMatrix& M);

  void GrammSchmidt(DenseMatrix& N);

 protected:
  int FindPosition_(int v, AmanziMesh::Entity_ID_List nodes);

  int stability_method_;  // stability parameters
  double scalar_stability_, scaling_factor_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

