/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The package uses the formula M = Mc + Ms, where matrix Mc is build from a
  consistency condition (Mc N = R) and matrix Ms is build from a stability
  condition and Ms N = 0, to generate mass and stiffness matrices for
  a variety of physics packages: flow, transport, thermal, and geomechanics.
  The material properties are imbedded into the the matrix Mc.

  Notation used below: M (mass), W (inverse of M), A (stiffness).
*/

#ifndef AMANZI_WHETSTONE_MFD3D_DIFFUSION_HH_
#define AMANZI_WHETSTONE_MFD3D_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "DeRham_Face.hh"
#include "Tensor.hh"
#include "MFD3D.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Diffusion : public DeRham_Face {
 public:
  // constructor for backward compatibility
  MFD3D_Diffusion(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh) : DeRham_Face(mesh){};
  MFD3D_Diffusion(const Teuchos::ParameterList& plist,
                  const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
    : DeRham_Face(mesh){};

  // main methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override
  {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, 1));
  }

  // -- default Derahm complex for the mass matrix is not used by Amanzi
  //    but we have to tell compiler a proper member function
  using DeRham_Face::MassMatrix;

  // -- inverse mass matrix is modified to reflect scaling of fluxes by area
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W) override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // -- divergence matrix
  virtual int DivergenceMatrix(int c, DenseMatrix& A) override;

  // other mimetic methods
  // -- bad consistency conditions (flux is scaled by area)
  int
  L2consistencyScaledArea(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  int L2consistencyInverseScaledArea(int c,
                                     const Tensor& K,
                                     DenseMatrix& R,
                                     DenseMatrix& Wc,
                                     bool symmetry);
  int MassMatrixScaledArea(int c, const Tensor& K, DenseMatrix& M);

  // -- optimized stability
  int MassMatrixInverseOptimized(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseMMatrixHex(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseMMatrix(int c, const Tensor& K, DenseMatrix& W);

  int StiffnessMatrixOptimized(int c, const Tensor& K, DenseMatrix& A);
  int StiffnessMatrixMMatrix(int c, const Tensor& K, DenseMatrix& A);

  // -- tensor is product k K
  int L2consistencyInverseDivKScaled(int c,
                                     const Tensor& K,
                                     double kmean,
                                     const AmanziGeometry::Point& kgrad,
                                     DenseMatrix& R,
                                     DenseMatrix& Wc);
  int MassMatrixInverseDivKScaled(int c,
                                  const Tensor& K,
                                  double kmean,
                                  const AmanziGeometry::Point& kgrad,
                                  DenseMatrix& W);

  // -- non-symmetric tensor K (consistency is not changed)
  int MassMatrixNonSymmetric(int c, const Tensor& K, DenseMatrix& M);
  int MassMatrixInverseNonSymmetric(int c, const Tensor& K, DenseMatrix& W);

  // surface methods
  // -- mass matrix
  int L2consistencyInverseSurface(int c, const Tensor& K, DenseMatrix& R, DenseMatrix& Wc);
  int MassMatrixInverseSurface(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseSurfaceTPFA(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseSurfaceMMatrix(int c, const Tensor& K, DenseMatrix& W);

  // -- other related discetization methods
  int MassMatrixInverseSO(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseTPFA(int c, const Tensor& K, DenseMatrix& W);
  int MassMatrixInverseDiagonal(int c, const Tensor& K, DenseMatrix& W);

  // -- projectors
  //    we return linear polynomial instead of constant vector polynomial (FIXME)
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& vc) override;

  // utils
  double Transmissibility(int f, int c, const Tensor& K);

 private:
  // stability methods (add stability matrix, M += Mstab)
  int StabilityMMatrixHex_(int c, const Tensor& K, DenseMatrix& M);
  void RescaleMassMatrixInverse_(int c, DenseMatrix& W);

  // mesh extension methods
  // -- exterior normal
  AmanziGeometry::Point mesh_face_normal(int f, int c);

 protected:
  using MFD3D::mesh_;
  using MFD3D::d_;

 private:
  static RegisteredFactory<MFD3D_Diffusion> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
