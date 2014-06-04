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
#include "Operator.hh"
#include "OperatorDiffusion.hh"


namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionWithGravity : public OperatorDiffusion {
 public:
  OperatorDiffusionWithGravity() {};
  OperatorDiffusionWithGravity(Teuchos::RCP<const CompositeVectorSpace> cvs, 
                               Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc) 
      : OperatorDiffusion(cvs, plist, bc) {};
  OperatorDiffusionWithGravity(const Operator& op, 
                               Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc) 
      : OperatorDiffusion(op, plist, bc) {};

  ~OperatorDiffusionWithGravity() {};

  // main members
  void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux, Teuchos::RCP<const CompositeVector> u);
  void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }

 private:
  AmanziGeometry::Point g_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

