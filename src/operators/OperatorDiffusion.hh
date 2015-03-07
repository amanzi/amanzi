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


/*
  Pure interface for Diffusion operators  
*/ 

namespace Amanzi {
namespace Operators {

class OperatorDiffusion {
 public:
  OperatorDiffusion(const Teuchos::RCP<Operator>& global_op) :
      global_op_(global_op) {};

  OperatorDiffusion(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      mesh_(mesh) {};

  OperatorDiffusion(const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      mesh_(mesh) {};
  
  // main virtual members
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     double rho, double mu) = 0;
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     const Teuchos::RCP<const CompositeVector>& rho,
                     const Teuchos::RCP<const CompositeVector>& mu) = 0;
  virtual void Setup(const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp) = 0;

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u) = 0;
  virtual void UpdateFlux(const CompositeVector& u, CompositeVector& flux) = 0;
  virtual void ApplyBCs(bool primary = true) = 0;
  virtual void ModifyMatrices(const CompositeVector& u) = 0;

  // default implementation  
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) {
    Setup(K, 1.0, 1.0);
  }
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp,
                     double rho, double mu) {
    Setup(K, rho, mu);
    Setup(k, dkdp);
  }
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp,
                     const Teuchos::RCP<const CompositeVector>& rho,
                     const Teuchos::RCP<const CompositeVector>& mu) {
    Setup(K, rho, mu);
    Setup(k, dkdp);
  }

  // boundary conditions
  virtual void SetBCs(const Teuchos::RCP<BCs>& bc) {
    bc_ = bc;
    global_op_->SetBCs(bc);
  }

  // gravity terms -- may not be implemented
  virtual void SetGravity(const AmanziGeometry::Point& g) {
    Errors::Message msg("OperatorDiffusion: This diffusion implementation does not support gravity.");
    Exceptions::amanzi_throw(msg);
  }
  virtual void SetGravityDensity(const Teuchos::RCP<const CompositeVector>& rho) {
    Errors::Message msg("OperatorDiffusion: This diffusion implementation does not support gravity.");
    Exceptions::amanzi_throw(msg);
  }
  
  // access
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  int schema_prec_dofs() { return global_op_schema_; }

  Teuchos::RCP<const Op> local_matrices() const { return local_op_; }
  Teuchos::RCP<Op> local_matrices() { return local_op_; }
  int schema_dofs() { return local_op_schema_; }

  Teuchos::RCP<const Op> jacobian_matrices() const { return jac_op_; }
  Teuchos::RCP<Op> jacobian_matrices() { return jac_op_; }
  int schema_jacobian() { return jac_op_schema_; }
  
  int upwind() { return upwind_; }
  
 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;

  // physics
  bool scalar_rho_mu_;
  double rho_, mu_;
  Teuchos::RCP<const CompositeVector> rho_cv_, mu_cv_;

  Teuchos::RCP<const CompositeVector> k_, dkdp_;

  int upwind_;

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


