/*
  This is the Operator component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete diffusion operator of a surface.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_SURFACE_HH_
#define AMANZI_OPERATOR_DIFFUSION_SURFACE_HH_

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

class OperatorDiffusionSurface : public Operator {
 public:
  OperatorDiffusionSurface() {};
  OperatorDiffusionSurface(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) : Operator(cvs, dummy) {};
  OperatorDiffusionSurface(const Operator& op) : Operator(op) {};
  ~OperatorDiffusionSurface() {};

  // main members
  void InitOperator(std::vector<WhetStone::Tensor>& K, Teuchos::RCP<NonlinearCoefficient> k,
                    const Teuchos::ParameterList& plist);
  void AssembleMatrix(int schema);
  void UpdateMatrices();

  // local implementation of matrix inversion
  int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;

  // local preconditioner
  void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist,
                          std::vector<int>& bc_model, std::vector<double>& bc_values);

 private:
  void CreateMassMatrices_(std::vector<WhetStone::Tensor>& K);

 private:
  Teuchos::ParameterList plist_;
  std::vector<WhetStone::DenseMatrix> Wff_cells_;

  int upwind_;
  Teuchos::RCP<NonlinearCoefficient> k_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


