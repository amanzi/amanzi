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

#ifndef AMANZI_WHETSTONE_MFD3D_DIFFUSION_HH_
#define AMANZI_WHETSTONE_MFD3D_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "DeRham_Face.hh"
#include "Tensor.hh"
#include "MFD3D.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Diffusion : public virtual MFD3D,
                        public virtual DeRham_Face { 
 public:
  explicit MFD3D_Diffusion(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh), 
      InnerProduct(mesh) {};
  ~MFD3D_Diffusion() {};

  // main methods 
  // -- default Derahm complex for the mass matrix is not used by Amanzi
  // -- inverse mass matrix is adjusted to reflect scaling of fluxes by area
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W); 

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A);

  // -- divergence matrix
  virtual int DivergenceMatrix(int c, DenseMatrix& A);

  // -- other interfaces not used yet by this class
  virtual int MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M) {};
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A) {};

  virtual int StiffnessMatrixPoly(int c, const Polynomial& K, DenseMatrix& A) {};
  virtual int AdvectionMatrixPoly(int c, const VectorPolynomial& v, DenseMatrix& A) {};

  // other mimetic methods
  // -- bad consistency conditions (flux is scaled by area)
  int L2consistencyScaledArea(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  int L2consistencyInverseScaledArea(int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);
  int MassMatrixScaledArea(int c, const Tensor& K, DenseMatrix& M);

  // -- optimized stability
  int MassMatrixInverseOptimized(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseMMatrixHex(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseMMatrix(int c, const Tensor& K, DenseMatrix& W);

  int StiffnessMatrixOptimized(int c, const Tensor& K, DenseMatrix& A);
  int StiffnessMatrixMMatrix(int c, const Tensor& K, DenseMatrix& A);

  // -- edge-based degrees of freedom
  int H1consistencyEdge(int cell, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int StiffnessMatrixEdge(int c, const Tensor& K, DenseMatrix& A);

  // -- tensor is product k K
  int L2consistencyInverseDivKScaled(int c, const Tensor& K,
                                     double kmean, const AmanziGeometry::Point& kgrad,
                                     DenseMatrix& R, DenseMatrix& Wc);
  int MassMatrixInverseDivKScaled(int c, const Tensor& K,
                                  double kmean, const AmanziGeometry::Point& kgrad, DenseMatrix& W);

  // -- non-symmetric tensor K (consistency is not changed)
  int MassMatrixNonSymmetric(int c, const Tensor& K, DenseMatrix& M);
  int MassMatrixInverseNonSymmetric(int c, const Tensor& K, DenseMatrix& W);

  // surface methods
  // -- mass matrix
  int L2consistencyInverseSurface(int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc);
  int MassMatrixInverseSurface(int c, const Tensor& K, DenseMatrix& W);

  // other related discetization methods
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

  // mesh extension methods 
  // -- exterior normal
  AmanziGeometry::Point mesh_face_normal(int f, int c);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

