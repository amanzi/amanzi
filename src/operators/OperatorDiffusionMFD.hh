/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_HH_
#define AMANZI_OPERATOR_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusion {
 public:
  OperatorDiffusion(Teuchos::ParameterList& plist,
                    Teuchos::RCP<Operator> global_op) :
      global_op_(global_op),
      mesh_(Teuchos::null),
      factor_(1.0)
  {
    InitDiffusion_(plist);
  }

  OperatorDiffusion(Teuchos::ParameterList& plist,
                    Teuchos::RCP<AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh),      
      factor_(1.0)
  {
    InitDiffusion_(plist);
  }

  OperatorDiffusion(Teuchos::ParameterList& plist,
                    Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh),      
      factor_(1.0)
  {
    InitDiffusion_(plist);
  }

  // main members
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K, double rho, double mu);
  virtual void Setup(Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp);
  virtual void Setup(Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                     Teuchos::RCP<const CompositeVector> rho, Teuchos::RCP<const CompositeVector> mu);

  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                     double rho, double mu);
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
                     Teuchos::RCP<const CompositeVector> rho, Teuchos::RCP<const CompositeVector> mu);

  virtual void UpdateMatrices(Teuchos::RCP<const CompositeVector> flux, Teuchos::RCP<const CompositeVector> u);
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux);
  virtual void ApplyBCs(const Teuchos::RCP<BCs>& bc, bool primary=true);

  // access (for developers mainly)
  void set_factor(double factor) { factor_ = factor; }
  int schema_dofs() { return local_op_schema_; }
  int schema_prec_dofs() { return global_op_schema_; }

  int nfailed_primary() { return nfailed_primary_; }

  // special members
  void ModifyMatrices(const CompositeVector& u);

  // access
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  Teuchos::RCP<const Op> local_matrices() const { return local_op_; }
  Teuchos::RCP<Op> local_matrices() { return local_op_; }

  int upwind() { return upwind_; }
  
 protected:
  void InitDiffusion_(Teuchos::ParameterList& plist);
  void CreateMassMatrices_();

  void UpdateMatricesNodal_();
  void UpdateMatricesTPFA_();
  void UpdateMatricesMixed_(Teuchos::RCP<const CompositeVector> flux);
  void UpdateMatricesMixedWithGrad_(Teuchos::RCP<const CompositeVector> flux);

  void AddNewtonCorrectionCell_(Teuchos::RCP<const CompositeVector> flux,
          Teuchos::RCP<const CompositeVector> u);

 protected:
  std::vector<WhetStone::DenseMatrix> Wff_cells_;
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  double rho_, mu_;
  bool scaled_constraint_;
  Teuchos::RCP<const CompositeVector> rho_cv_, mu_cv_;

  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  int upwind_;
  int newton_correction_;
  double factor_;

  int mfd_primary_, mfd_secondary_, mfd_pc_primary_, mfd_pc_secondary_;
  int nfailed_primary_;
  bool scalar_rho_mu_;

  // operator
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  Teuchos::RCP<Op> jac_op_;
  int global_op_schema_, local_op_schema_, jac_op_schema_;
  Teuchos::RCP<BCs> bc_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;

};

}  // namespace Operators
}  // namespace Amanzi


#endif


