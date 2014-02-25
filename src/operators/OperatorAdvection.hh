/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete advection operator of a surface.
*/

#ifndef AMANZI_OPERATOR_ADVECTION_HH_
#define AMANZI_OPERATOR_ADVECTION_HH_

#include "Epetra_IntVector.h"

#include "Operator.hh"


namespace Amanzi {
namespace Operators {

class OperatorAdvection : public Operator {
 public:
  OperatorAdvection() {};
  OperatorAdvection(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorAdvection(const Operator& op) : Operator(op) {};
  ~OperatorAdvection() {};

  // main members
  void InitOperator(const CompositeVector& u);
  void UpdateMatrices(const CompositeVector& u);

 private:
  void IdentifyUpwindCells_(const CompositeVector& u);

 private:
  Teuchos::RCP<Epetra_IntVector> upwind_cell_, downwind_cell_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

