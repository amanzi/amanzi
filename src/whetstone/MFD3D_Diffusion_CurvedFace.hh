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

#ifndef AMANZI_WHETSTONE_MFD3D_DIFFUSION_CURVED_FACE_HH_
#define AMANZI_WHETSTONE_MFD3D_DIFFUSION_CURVED_FACE_HH_

#include <memory>

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

class MFD3D_Diffusion_CurvedFace : public DeRham_Face {
 public:
  // constructor for backward compatibility
  MFD3D_Diffusion_CurvedFace(const Teuchos::ParameterList& plist,
                             const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
    : DeRham_Face(mesh){};
  ~MFD3D_Diffusion_CurvedFace(){};

  // main methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override
  {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, d_));
  }

  // -- mass matrix
  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override;

  // -- inverse mass matrix is modified to reflect scaling of fluxes by area
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W) override;

  // -- stiffness matrix
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // -- divergence matrix
  virtual int DivergenceMatrix(int c, DenseMatrix& A) override { return 0; };

  // -- projectors
  //    we return linear polynomial instead of constant vector polynomial (FIXME)
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& vc) override;

  void set_generalized_centroids(const std::shared_ptr<std::vector<AmanziGeometry::Point>>& bf) { bf_ = bf; }

 private:
  void RescaleMassMatrixInverse_(int c, DenseMatrix& W);

 protected:
  using MFD3D::mesh_;
  using MFD3D::d_;

  std::shared_ptr<const std::vector<AmanziGeometry::Point>> bf_;

 private:
  static RegisteredFactory<MFD3D_Diffusion_CurvedFace> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
