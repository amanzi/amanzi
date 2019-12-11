/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! Diffusion generates local Ops and global Operators for an elliptic operator.

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
#include "PDE_HelperDiscretization.hh"

/*!
Example:

.. code-block:: xml

    <ParameterList name="OPERATOR_NAME">
      <Parameter name="discretization primary" type="string"
                value="mfd: optimized for monotonicity"/>
      <Parameter name="discretization secondary" type="string"
                value="mfd: two-point flux approximation"/>
      <Parameter name="schema" type="Array(string)"
                value="{face, cell}"/>
      <Parameter name="preconditioner schema" type="Array(string)"
                value="{face}"/>
      <Parameter name="gravity" type="bool" value="true"/>
      <Parameter name="gravity term discretization" type="string"
                value="hydraulic head"/>
      <Parameter name="nonlinear coefficient" type="string"
                value="upwind: face"/>
      <Parameter name="Newton correction" type="string"
                value="true Jacobian"/>

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


/*
  Ghost elemets of composite vectors k_ and dkdp_ are NOT used to be up to date.
  They always have to be updated from master elements before accessing their
  values.
*/

namespace Amanzi {
namespace Operators {

class PDE_Diffusion : public PDE_HelperDiscretization {
 public:
  // using TensorCoef_view_type = Kokkos::View<double**>;
  using TensorCoef_view_type = Kokkos::View<double**, Kokkos::LayoutRight>;
  
  PDE_Diffusion(const Teuchos::RCP<Operator>& global_op)
    : PDE_HelperDiscretization(global_op),
      K_(),
      k_(Teuchos::null),
      dkdp_(Teuchos::null){};

  PDE_Diffusion(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh),
      K_(),
      k_(Teuchos::null),
      dkdp_(Teuchos::null){};

  PDE_Diffusion(const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
    : PDE_HelperDiscretization(mesh),
      K_(),
      k_(Teuchos::null),
      dkdp_(Teuchos::null){};

  virtual ~PDE_Diffusion() = default;

  // main virtual members
  // -- setup
  virtual void SetTensorCoefficient(
      const WhetStone::TensorArray& K) = 0;
  virtual void
  SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                       const Teuchos::RCP<const CompositeVector>& dkdp) = 0;

  // -- creation of an operator
  virtual void UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    double scalar_factor = 1.0) = 0;

  virtual void UpdateMatricesNewtonCorrection(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& factor) = 0;

  // -- additional interface on non-manifolds
  virtual void
  UpdateFluxNonManifold(const Teuchos::Ptr<const CompositeVector>& u,
                        const Teuchos::Ptr<CompositeVector>& flux) = 0;

  // -- matrix modifications
  virtual void ModifyMatrices(const CompositeVector& u) = 0;
  virtual void ScaleMassMatrices(double s) = 0;

  // default implementation
  virtual void
  Setup(const WhetStone::TensorArray& K,
        const Teuchos::RCP<const CompositeVector>& k,
        const Teuchos::RCP<const CompositeVector>& dkdp)
  {
    SetTensorCoefficient(K);
    SetScalarCoefficient(k, dkdp);
  }

  // -- working with consistent faces -- may not be implemented
  virtual int UpdateConsistentFaces(CompositeVector& u)
  {
    Errors::Message msg("Diffusion: This diffusion implementation does not "
                        "support working with consistent faces.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // interface to solvers for treating nonlinear BCs.
  virtual double ComputeTransmissibility(int f) const = 0;
  // virtual double ComputeGravityFlux(int f) const = 0;

  // access
  int schema_prec_dofs() { return global_op_schema_; }
  int schema_dofs() { return local_op_schema_; }

  Teuchos::RCP<const Op> jacobian_matrices() const { return jac_op_; }
  Teuchos::RCP<Op> jacobian_matrices() { return jac_op_; }
  void set_jacobian_matrices(const Teuchos::RCP<Op>& op)
  {
    if (global_operator().get()) {
      if (local_matrices().get()) {
        auto index = std::find(global_operator()->OpBegin(),
                               global_operator()->OpEnd(),
                               jac_op_) -
                     global_operator()->OpBegin();
        if (index != global_operator()->OpSize()) {
          global_operator()->OpPushBack(op);
        } else {
          global_operator()->OpReplace(op, index);
        }
      } else {
        global_operator()->OpPushBack(op);
      }
    }
    jac_op_ = op;
  }
  int schema_jacobian() { return jac_op_schema_; }

  int little_k() const { return little_k_; }
  CompositeVectorSpace little_k_space() const
  {
    CompositeVectorSpace out;
    out.SetMesh(mesh_);
    out.SetGhosted();
    if (little_k_ == OPERATOR_LITTLE_K_NONE) { return out; }
    if (little_k_ != OPERATOR_LITTLE_K_UPWIND) {
      out.AddComponent("cell", AmanziMesh::CELL, 1);
    }
    if (little_k_ != OPERATOR_LITTLE_K_STANDARD) {
      out.AddComponent("face", AmanziMesh::FACE, 1);
    }
    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN ||
        little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      out.AddComponent("twin", AmanziMesh::FACE, 1);
    }
    if (little_k_ == OPERATOR_LITTLE_K_DIVK_TWIN_GRAD) {
      out.AddComponent("grad", AmanziMesh::CELL, mesh_->space_dimension());
    }
    return out;
  }

 protected:
  WhetStone::TensorArray K_;
  bool K_symmetric_;

  // nonlinear coefficient and its representation
  Teuchos::RCP<const CompositeVector> k_, dkdp_;
  int little_k_;

  // additional operators
  Teuchos::RCP<Op> jac_op_;
  int global_op_schema_, local_op_schema_, jac_op_schema_;
};

} // namespace Operators
} // namespace Amanzi

#endif
