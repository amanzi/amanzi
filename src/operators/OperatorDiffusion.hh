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
  OperatorDiffusion() { InitDiffusion_(); };
  OperatorDiffusion(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) { InitDiffusion_(); }
  OperatorDiffusion(const Operator& op) : Operator(op) { InitDiffusion_(); }
  ~OperatorDiffusion() {};

  // main members
  void InitOperator(std::vector<WhetStone::Tensor>& K, Teuchos::RCP<NonlinearCoefficient> k,
                    int schema_base, int schema_dofs, const Teuchos::ParameterList& plist);
  void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux);
  void UpdateMatricesStiffness(std::vector<WhetStone::Tensor>& K);
  void UpdateFlux(const CompositeVector& u, CompositeVector& flux, double scalar);

  // re-implementation of basic operator virtual members
  int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;

  // special assembling routines
  void AssembleMatrixSpecial();
  int ApplyInverseSpecial(const CompositeVector& X, CompositeVector& Y) const;
  void InitPreconditionerSpecial(const std::string& prec_name, const Teuchos::ParameterList& plist,
                                 std::vector<int>& bc_model, std::vector<double>& bc_values);

  // internal memebers
  void CreateMassMatrices_(std::vector<WhetStone::Tensor>& K);

  // access (for developers only)
  void set_factor(double factor) { factor_ = factor; }

 private:
  void InitDiffusion_();

 public:
  Teuchos::ParameterList plist_;
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  Teuchos::RCP<NonlinearCoefficient> k_;

  int schema_base_, schema_dofs_, schema_;
  int upwind_;

  double factor_;

 private:
  mutable bool special_assembling_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


