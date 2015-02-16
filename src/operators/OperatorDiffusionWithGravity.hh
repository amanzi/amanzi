/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete gravity operator blended with the diffusion operator.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_DIFFUSION_WITH_GRAVITY_HH_

#include "Epetra_IntVector.h"

#include "tensor.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionWithGravity : public OperatorDiffusion {
 public:
  OperatorDiffusionWithGravity(Teuchos::ParameterList& plist,
                    Teuchos::RCP<Operator> global_op) :
      OperatorDiffusion(plist, global_op)
  {
    Init_();
  }

  OperatorDiffusionWithGravity(Teuchos::ParameterList& plist,
                    Teuchos::RCP<AmanziMesh::Mesh> mesh) :
      OperatorDiffusion(plist, mesh)
  {
    Init_();
  }

  OperatorDiffusionWithGravity(Teuchos::ParameterList& plist,
                    Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      OperatorDiffusion(plist, mesh)
  {
    Init_();
  }

  // main members
  virtual void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux, Teuchos::RCP<const CompositeVector> u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }
  

 protected:
  inline AmanziGeometry::Point GravitySpecialDirection_(int f) const;
  void Init_() { gravity_special_projection_ = (mfd_primary_ == WhetStone::DIFFUSION_TPFA); }

 protected:
  AmanziGeometry::Point g_;
  bool gravity_special_projection_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

