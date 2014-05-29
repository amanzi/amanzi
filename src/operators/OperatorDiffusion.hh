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
#include "Point.hh"
#include "CompositeVector.hh"

#include "Operator.hh"
#include "OperatorTypeDefs.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusion : public Operator {
 public:
  OperatorDiffusion() {};
  OperatorDiffusion(Teuchos::RCP<const CompositeVectorSpace> cvs, const Teuchos::ParameterList& plist) 
      : Operator(cvs, 0) { InitDiffusion_(plist); }
  OperatorDiffusion(const Operator& op, const Teuchos::ParameterList& plist)
      : Operator(op) { InitDiffusion_(plist); }
  ~OperatorDiffusion() {};

  // main members
  void InitOperator(std::vector<WhetStone::Tensor>& K,
                    Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                    double rho, double mu);
  void InitOperator(std::vector<WhetStone::Tensor>& K,
                    Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                    Teuchos::RCP<const CompositeVector> rho, Teuchos::RCP<const CompositeVector> mu);
  virtual void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  // re-implementation of basic operator virtual members
  void AssembleMatrix(int schema);
  int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;

  // access (for developers only)
  void set_factor(double factor) { factor_ = factor; }
  int schema_dofs() { return schema_dofs_; }
  int schema_prec_dofs() { return schema_prec_dofs_; }

  // preconditioners
  void SetPreconditionerBCs(std::vector<int>& bc_model, std::vector<double>& bc_values);
  void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist,
                          std::vector<int>& bc_model, std::vector<double>& bc_values);

  // special members
  void ModifyMatrix(const CompositeVector& u);

 protected:
  void CreateMassMatrices_();

  void InitDiffusion_(const Teuchos::ParameterList& plist);
  void UpdateMatricesNodal_();
  void UpdateMatricesTPFA_();
  void UpdateMatricesMixed_(Teuchos::RCP<const CompositeVector> flux);
  int ApplyInverseSpecial_(const CompositeVector& X, CompositeVector& Y) const;
  void InitPreconditionerSpecialFE_(const std::string& prec_name, const Teuchos::ParameterList& plist,
                                    std::vector<int>& bc_model, std::vector<double>& bc_values);
  void InitPreconditionerSpecialCRS_(const std::string& prec_name, const Teuchos::ParameterList& plist,
                                     std::vector<int>& bc_model, std::vector<double>& bc_values);

 public:
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  std::vector<WhetStone::Tensor>* K_;
  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  double rho_, mu_;
  Teuchos::RCP<const CompositeVector> rho_cv_, mu_cv_;

  int schema_base_, schema_dofs_, schema_;
  int schema_prec_dofs_;
  mutable bool special_assembling_;

  double factor_;

  int mfd_primary_, mfd_secondary_;
  bool scalar_rho_mu_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


