/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method.
*/

#ifndef AMANZI_MFD3D_HH_
#define AMANZI_MFD3D_HH_

/*
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

#include "WhetStoneDefs.hh"
#include "WhetStone_typedefs.hh"
#include "DenseMatrix.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D { 
 public:
  explicit MFD3D(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D() {};

  virtual int L2consistency(int cell, const Tensor& T,
                            DenseMatrix& N, DenseMatrix& Mc, bool symmetry) = 0;

  virtual int L2consistencyInverse(int cell, const Tensor& T,
                                   DenseMatrix& R, DenseMatrix& Wc, bool symmetry) = 0;

  virtual int H1consistency(int cell, const Tensor& T,
                            DenseMatrix& N, DenseMatrix& Mc) = 0;

  virtual int MassMatrix(int cell, const Tensor& T, DenseMatrix& M) = 0; 

  virtual int MassMatrixInverse(int cell, const Tensor& T, DenseMatrix& W) = 0; 

  virtual int StiffnessMatrix(int cell, const Tensor& T, DenseMatrix& A) = 0; 

  // experimental methods (for stability region analysis; unit test)
  double CalculateStabilityScalar(DenseMatrix& Mc);
  void ModifyStabilityScalingFactor(double factor);

  // geometry methods
  void PolygonCentroidWeights(
      const Entity_ID_List& nodes,
      double area, std::vector<double>& weights) const;

  // access members
  double scaling_factor() { return scaling_factor_; }
  double scalar_stability() { return scalar_stability_; }
  double simplex_functional() { return simplex_functional_; }
  int simplex_num_itrs() { return simplex_num_itrs_; }

  // extension of the mesh API (must be removed lipnikov@lanl.gov)
  int cell_get_face_adj_cell(int cell, int face);

 protected:
  // supporting stability methods (add matrix M += Mstab)
  // use R, Wc, W for the inverse matrix
  void StabilityScalar(int c, DenseMatrix& N, DenseMatrix& M);

  int StabilityOptimized(const Tensor& T, DenseMatrix& N, DenseMatrix& M);

  int StabilityMonotoneHex(int c, const Tensor& T, DenseMatrix& Mc, DenseMatrix& M);

  int StabilityMMatrix_(int c, DenseMatrix& N, DenseMatrix& M, 
                        int objective = WHETSTONE_SIMPLEX_FUNCTIONAL_SUMALL);

  int SimplexFindFeasibleSolution_(DenseMatrix& T, int m1, int m2, int m3, int* izrow, int* iypos);
  void SimplexPivotElement_(DenseMatrix& T, int kp, int* ip);
  void SimplexExchangeVariables_(DenseMatrix& T, int kp, int ip);

  void GrammSchmidt(DenseMatrix& N);

 protected:
  int FindPosition_(int v, Entity_ID_List nodes);

  int stability_method_;  // stability parameters
  double scalar_stability_, scaling_factor_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  double simplex_functional_;
  int simplex_num_itrs_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

