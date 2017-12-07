// Diffusion generates local Ops and global Operators for an elliptic operator.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "Point.hh"
#include "CompositeVector.hh"
#include "DenseMatrix.hh"

#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_Helper.hh"

/*!
Example:

.. code-block:: xml

    <ParameterList name="OPERATOR_NAME">
      <Parameter name="discretization primary" type="string" value="mfd: optimized for monotonicity"/>
      <Parameter name="discretization secondary" type="string" value="mfd: two-point flux approximation"/>
      <Parameter name="schema" type="Array(string)" value="{face, cell}"/>
      <Parameter name="preconditioner schema" type="Array(string)" value="{face}"/>
      <Parameter name="gravity" type="bool" value="true"/>
      <Parameter name="gravity term discretization" type="string" value="hydraulic head"/>
      <Parameter name="nonlinear coefficient" type="string" value="upwind: face"/>
      <Parameter name="Newton correction" type="string" value="true Jacobian"/>

      <ParameterList name="consistent faces">
        <ParameterList name="linear solver">
          ...
        </ParameterList>
        <ParameterList name="preconditioner">
          ...
        </ParameterList>
      </ParameterList>
    </ParameterList>
*/


namespace Amanzi {
namespace Operators {

class PDE_Diffusion : public PDE_Helper {
 public:
  PDE_Diffusion(const Teuchos::RCP<Operator>& global_op) :
      PDE_Helper(global_op),
      K_(Teuchos::null),
      k_(Teuchos::null),
      dkdp_(Teuchos::null)
  {};

  PDE_Diffusion(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Helper(mesh),
      K_(Teuchos::null),
      k_(Teuchos::null),
      dkdp_(Teuchos::null)
  {};

  PDE_Diffusion(const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
      PDE_Helper(mesh),
      K_(Teuchos::null),
      k_(Teuchos::null),
      dkdp_(Teuchos::null)
  {};

  virtual ~PDE_Diffusion() = default;
  
  // main virtual members
  // -- setup 
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) = 0;
  virtual void SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                    const Teuchos::RCP<const CompositeVector>& dkdp) = 0;

  // -- creation of an operator
  virtual void UpdateMatricesNewtonCorrection(
          const Teuchos::Ptr<const CompositeVector>& flux,
          const Teuchos::Ptr<const CompositeVector>& u,
          double scalar_limiter = 1.0) = 0;

  // -- additional interface on non-manifolds
  virtual void UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
                                     const Teuchos::Ptr<CompositeVector>& flux) = 0;

  // -- matrix modification
  virtual void ApplyBCs(bool primary, bool eliminate) = 0;
  virtual void ModifyMatrices(const CompositeVector& u) = 0;
  virtual void ScaleMassMatrices(double s) = 0;

  // default implementation  
  virtual void Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K,
                     const Teuchos::RCP<const CompositeVector>& k,
                     const Teuchos::RCP<const CompositeVector>& dkdp) {
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
  }

  // boundary conditions (BC) require information on test and
  // trial spaces. For a single PDE, these BCs could be the same.
  virtual void SetBCs(const Teuchos::RCP<BCs>& bc_trial,
                      const Teuchos::RCP<BCs>& bc_test) {
    SetTrialBCs(bc_trial);  
    SetTestBCs(bc_test);  
  }
  virtual void SetTrialBCs(const Teuchos::RCP<BCs>& bc) {
    if (bcs_trial_.size() == 0) {
      bcs_trial_.resize(1);
    }
    bcs_trial_[0] = bc;
    global_op_->SetTrialBCs(bc);
  }
  virtual void SetTestBCs(const Teuchos::RCP<BCs>& bc) {
    if (bcs_test_.size() == 0) {
      bcs_test_.resize(1);
    }
    bcs_test_[0] = bc;
    global_op_->SetTestBCs(bc);
  }

  // -- working with consistent faces -- may not be implemented
  virtual int UpdateConsistentFaces(CompositeVector& u) {
    Errors::Message msg("Diffusion: This diffusion implementation does not support working with consistent faces.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
  
  // interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const = 0;
  virtual double ComputeGravityFlux(int f) const = 0;

  // access
  int schema_prec_dofs() { return global_op_schema_; }

  Teuchos::RCP<const Op> local_matrices() const { return local_op_; }
  Teuchos::RCP<Op> local_matrices() { return local_op_; }
  int schema_dofs() { return local_op_schema_; }

  Teuchos::RCP<const Op> jacobian_matrices() const { return jac_op_; }
  Teuchos::RCP<Op> jacobian_matrices() { return jac_op_; }
  int schema_jacobian() { return jac_op_schema_; }

  int little_k() { return little_k_; }
  
 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  bool K_symmetric_;

  // nonlinear coefficient and its representation
  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  int little_k_;

  // additional operators
  Teuchos::RCP<Op> jac_op_;
  int global_op_schema_, local_op_schema_, jac_op_schema_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


