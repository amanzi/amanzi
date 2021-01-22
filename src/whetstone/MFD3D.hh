/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method.

  The package uses the formula M = Mc + Ms, where matrix Mc is build from a 
  consistency condition (Mc N = R) and matrix Ms is build from a stability 
  condition (Ms N = 0), to generate mass and stiffness matrices for a variety 
  of physics packages: flow, transport, thermal, and geomechanics. 
  The material properties are imbedded into the the matrix Mc. 

  Notation used below: M (mass), W (inverse of M), A (stiffness).

  NOTE: This class should be never instantiated directly. It is used to
  add additional functionality specific for MFD methods, such as various
  stabilization terms.
*/

#ifndef AMANZI_MFD3D_HELPER_HH_
#define AMANZI_MFD3D_HELPER_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "BilinearForm.hh"
#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D : public BilinearForm {
 public:
  MFD3D(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);

  // access members
  double simplex_functional() { return simplex_functional_; }
  int simplex_num_itrs() { return simplex_num_itrs_; }

  // experimental methods (for stability region analysis; unit test)
  void ModifyStabilityScalingFactor(double factor);

 protected:
  void StabilityScalar_(DenseMatrix& N, DenseMatrix& M);
  void StabilityScalarNonSymmetric_(DenseMatrix& N, DenseMatrix& M);
  int StabilityOptimized_(const Tensor& T, DenseMatrix& N, DenseMatrix& M);

  double CalculateStabilityScalar_(DenseMatrix& Mc);
  void GrammSchmidt_(DenseMatrix& N);

  // supporting stability methods (add matrix M += Mstab)
  // use R, Wc, W for the inverse matrix
  int StabilityMonotoneHex(int c, const Tensor& T, DenseMatrix& Mc, DenseMatrix& M);

  int StabilityMMatrix_(int c, DenseMatrix& N, DenseMatrix& M, 
                        int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_SUMALL);

  int SimplexFindFeasibleSolution_(DenseMatrix& T, int m1, int m2, int m3, int* izrow, int* iypos);
  void SimplexPivotElement_(DenseMatrix& T, int kp, int* ip);
  void SimplexExchangeVariables_(DenseMatrix& T, int kp, int ip);

 protected:
  int stability_method_;  // stability parameters
  double scalar_stability_, scaling_factor_;

  double simplex_functional_;
  int simplex_num_itrs_;
};


// non-member functions
void AddGradient(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh, int c, DenseMatrix& N);

}  // namespace WhetStone
}  // namespace Amanzi

#endif

