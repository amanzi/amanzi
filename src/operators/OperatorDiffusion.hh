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
#include "OperatorDefs.hh"


namespace Amanzi {
namespace Operators {

class OperatorDiffusion : public Operator {
 public:
  OperatorDiffusion() {};
  OperatorDiffusion(Teuchos::RCP<const CompositeVectorSpace> cvs, 
                    Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc) 
      : Operator(cvs, 0) { InitDiffusion_(bc, plist); }
  OperatorDiffusion(const Operator& op, Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc)
      : Operator(op) { InitDiffusion_(bc, plist); }
  ~OperatorDiffusion() {};

  // main members
  virtual void InitOperator(std::vector<WhetStone::Tensor>& K,
                            Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                            double rho, double mu);
  virtual void InitOperator(std::vector<WhetStone::Tensor>& K,
                            Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                            Teuchos::RCP<const CompositeVector> rho, Teuchos::RCP<const CompositeVector> mu);

  virtual void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux, Teuchos::RCP<const CompositeVector> u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);

  template <class WRM> 
  double DeriveBoundaryFaceValue(int f, const CompositeVector& u, const WRM&){

    if (u.HasComponent("face")) {
      const Epetra_MultiVector& u_face = *u.ViewComponent("face");
      return u_face[f][0];
    } else {
      const std::vector<int>& bc_model = bc_->bc_model();
      const std::vector<double>& bc_value = bc_->bc_value();
      if (bc_model[f] == OPERATOR_BC_FACE_DIRICHLET){
	return bc_value[f];
      } else {
	const Epetra_MultiVector& u_cell = *u.ViewComponent("cell");
	AmanziMesh::Entity_ID_List cells;
	mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
	int c = cells[0];
	return u_cell[0][c];
      }
    }
  };



  // re-implementation of basic operator virtual members
  void AssembleMatrix(int schema);
  int ApplyInverse(const CompositeVector& X, CompositeVector& Y) const;

  void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist);

  // access (for developers only)
  void set_factor(double factor) { factor_ = factor; }
  int schema_dofs() { return schema_dofs_; }
  int schema_prec_dofs() { return schema_prec_dofs_; }

  // special members
  void ModifyMatrices(const CompositeVector& u);

  // access
  int nfailed_primary() { return nfailed_primary_; }

 protected:
  void CreateMassMatrices_();

  void InitDiffusion_(Teuchos::RCP<BCs> bc, Teuchos::ParameterList& plist);
  void UpdateMatricesNodal_();
  void UpdateMatricesTPFA_();
  void UpdateMatricesMixed_(Teuchos::RCP<const CompositeVector> flux);
  int ApplyInverseSpecialSff_(const CompositeVector& X, CompositeVector& Y) const;
  int ApplyInverseSpecialScc_(const CompositeVector& X, CompositeVector& Y) const;
  void InitPreconditionerSpecialSff_(const std::string& prec_name, const Teuchos::ParameterList& plist);
  void InitPreconditionerSpecialScc_(const std::string& prec_name, const Teuchos::ParameterList& plist);

 public:
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  std::vector<WhetStone::Tensor>* K_;
  double rho_, mu_;
  Teuchos::RCP<const CompositeVector> rho_cv_, mu_cv_;

  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  int upwind_;

  int schema_base_, schema_dofs_, schema_;
  int schema_prec_dofs_;
  mutable int special_assembling_;

  double factor_;

  int mfd_primary_, mfd_secondary_;
  int nfailed_primary_;
  bool scalar_rho_mu_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


