/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete accumulation operator.
*/

#ifndef AMANZI_OPERATOR_ACCUMULATION_HH_
#define AMANZI_OPERATOR_ACCUMULATION_HH_

#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class OperatorAccumulation : public Operator {
 public:
  OperatorAccumulation() {};
  OperatorAccumulation(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorAccumulation(const Operator& op) : Operator(op) {};
  ~OperatorAccumulation() {};

  // main members
  void UpdateMatrices(const CompositeVector& u0, const CompositeVector& ss, double dT);
};

}  // namespace Operators
}  // namespace Amanzi

#endif

