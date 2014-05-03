/*
  This is the Operator component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete diffusion operator.
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
#include "NonlinearCoefficient.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusion : public Operator {
 public:
  OperatorDiffusion() : factor_(1.0) {};
  OperatorDiffusion(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorDiffusion(const Operator& op) : Operator(op) {};
  ~OperatorDiffusion() {};

  // main members
  void InitOperator(std::vector<WhetStone::Tensor>& K, Teuchos::RCP<NonlinearCoefficient> k,
                    int schema_base, int schema_dofs, const Teuchos::ParameterList& plist);
  void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux);
  void UpdateMatricesStiffness(std::vector<WhetStone::Tensor>& K);
  void UpdateFlux(const CompositeVector& u, CompositeVector& flux, double scalar);

  // access (for developers only)
  void set_factor(double factor) { factor_ = factor; }

 private:
  void CreateMassMatrices_(std::vector<WhetStone::Tensor>& K);

 private:
  Teuchos::ParameterList plist_;
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  Teuchos::RCP<NonlinearCoefficient> k_;

  int schema_base_, schema_dofs_, schema_;
  double factor_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


