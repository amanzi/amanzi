/*
  This is the Operator component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete diffusion operator of a surface.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_HH_
#define AMANZI_OPERATOR_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "tensor.hh"
#include "CompositeVector.hh"

#include "Operator.hh"
#include "OperatorTypeDefs.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusion : public Operator {
 public:
  OperatorDiffusion() {};
  OperatorDiffusion(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorDiffusion(const Operator& op) : Operator(op) {};
  ~OperatorDiffusion() {};

  // main members
  void UpdateMatricesStiffness(std::vector<WhetStone::Tensor>& K);
};

}  // namespace Operators
}  // namespace Amanzi


#endif


