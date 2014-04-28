/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete gravity operator.
*/

#ifndef AMANZI_OPERATOR_GRAVITY_HH_
#define AMANZI_OPERATOR_GRAVITY_HH_

#include "Epetra_IntVector.h"

#include "tensor.hh"
#include "Operator.hh"
#include "NonlinearCoefficient.hh"


namespace Amanzi {
namespace Operators {

class OperatorGravity : public Operator {
 public:
  OperatorGravity() {};
  OperatorGravity(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorGravity(const Operator& op) : Operator(op) {};
  ~OperatorGravity() {};

  // main members
  void InitOperator(Teuchos::RCP<NonlinearCoefficient> k) { k_ = k; }
  void UpdateMatrices(const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& rho_g);
  void UpdateFlux(const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& rho_g, 
                  CompositeVector& flux, double scalar);

 private:
  Teuchos::RCP<NonlinearCoefficient> k_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

