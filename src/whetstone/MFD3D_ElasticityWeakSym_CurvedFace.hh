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

  The mimetic finite difference method for elasticity with weak symmetry.
  Stress: 1 normal component per face (d dofs)
  Displacement: 1 value per cell (d dofs)
  Rotations:    1 value per cell (d dofs in 3D, 1 dof in 2D)
*/

#ifndef AMANZI_MFD3D_ELASTICITY_WEAK_SYMMETRY_CURVED_FACE_HH_
#define AMANZI_MFD3D_ELASTICITY_WEAK_SYMMETRY_CURVED_FACE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "MFD3D_ElasticityWeakSym.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_ElasticityWeakSym_CurvedFace : public MFD3D_ElasticityWeakSym {
 public:
  MFD3D_ElasticityWeakSym_CurvedFace(const Teuchos::ParameterList& plist,
                                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D_ElasticityWeakSym(plist, mesh) {};

  // mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc) override;

  // rotation metrix
  virtual void RotationMatrix(int c, DenseMatrix& G) override;

  void set_generalized_centroids(const std::shared_ptr<std::vector<AmanziGeometry::Point>>& bf)
  {
    bf_ = bf;
  }

 private:
  std::shared_ptr<const std::vector<AmanziGeometry::Point>> bf_;

  static RegisteredFactory<MFD3D_ElasticityWeakSym_CurvedFace> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
